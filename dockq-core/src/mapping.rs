//! Chain grouping, mapping enumeration, and the parallel multimer search.
//! Ports `format_mapping`, `group_chains`, `product_without_dupl` (exact enumeration
//! order), `get_all_chain_maps`, `run_on_all_native_interfaces`, and `main`'s search.
//!
//! Performance: every distinct interface (model-pair × native-pair) is scored exactly
//! once, in parallel (Rayon), then mappings are assembled from that cache — mirroring the
//! reference's `@lru_cache` on `run_on_chains` but shared across the whole search. The
//! best mapping is chosen by a deterministic argmax (max total DockQ, first in
//! enumeration order on ties — matching the reference's strict `>`), so results are
//! independent of thread scheduling.

use std::collections::HashMap;

use indexmap::{IndexMap, IndexSet};
use rayon::prelude::*;

use crate::align;
use crate::dockq;
use crate::error::{DockQError, Result};
use crate::model::{Chain, InterfaceResult, Structure};

/// Options for a single model-vs-native scoring run (mirrors the CLI flags).
#[derive(Clone, Debug)]
pub struct RunOptions {
    pub no_align: bool,
    pub capri_peptide: bool,
    pub small_molecule: bool,
    /// `--mapping MODEL:NATIVE` string (with optional `*` wildcards), if any.
    pub mapping: Option<String>,
    pub allowed_mismatches: usize,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self {
            no_align: false,
            capri_peptide: false,
            small_molecule: false,
            mapping: None,
            allowed_mismatches: 0,
        }
    }
}

/// Full result of a single model-vs-native run (mirrors `main`'s `info`, minus printing).
#[derive(Clone, Debug)]
pub struct RunResult {
    pub model: String,
    pub native: String,
    /// `models:natives` string in best-mapping order.
    pub best_mapping_str: String,
    /// Sum of per-interface DockQ over the best mapping.
    pub best_dockq: f64,
    /// `best_dockq / number_of_interfaces`.
    pub global_dockq: f64,
    /// Per-interface results, keyed by native chain pair (e.g. "AB"), in scoring order.
    pub best_result: IndexMap<String, InterfaceResult>,
    /// Best native→model chain map.
    pub best_mapping: IndexMap<String, String>,
}

/// Align the two chain pairs and run the per-interface calculation (ports `run_on_chains`).
/// Small-molecule interfaces hard-error (deferred; no silent fallback).
fn run_on_chains(
    model_chains: (&Chain, &Chain),
    native_chains: (&Chain, &Chain),
    no_align: bool,
    capri_peptide: bool,
    small_molecule: bool,
) -> Result<Option<InterfaceResult>> {
    if small_molecule {
        return Err(DockQError::SmallMoleculeUnsupported);
    }
    let aln0 = align::align_chains(model_chains.0, native_chains.0, no_align);
    let aln1 = align::align_chains(model_chains.1, native_chains.1, no_align);
    dockq::calc_dockq(model_chains, native_chains, (&aln0, &aln1), capri_peptide)
}

/// Score every native interface for an explicit native→model chain map (ports the public
/// `run_on_all_native_interfaces`). This is the single-map path used by the drop-in Python
/// API; the multimer search uses a cached variant internally. Returns (per-interface
/// results keyed by native chain pair, total DockQ).
pub fn run_on_native_interfaces(
    model: &Structure,
    native: &Structure,
    chain_map: &IndexMap<String, String>,
    no_align: bool,
    capri_peptide: bool,
) -> Result<(IndexMap<String, InterfaceResult>, f64)> {
    let mut result_mapping: IndexMap<String, InterfaceResult> = IndexMap::new();
    let native_ids: Vec<&String> = chain_map.keys().collect();

    for a in 0..native_ids.len() {
        for b in (a + 1)..native_ids.len() {
            let (na, nb) = (native_ids[a], native_ids[b]);
            let (ma, mb) = (&chain_map[na], &chain_map[nb]);
            if ma == mb {
                continue;
            }
            let nc0 = native
                .chain(na)
                .ok_or_else(|| DockQError::ChainNotFound(na.clone()))?;
            let nc1 = native
                .chain(nb)
                .ok_or_else(|| DockQError::ChainNotFound(nb.clone()))?;
            let mc0 = model
                .chain(ma)
                .ok_or_else(|| DockQError::ChainNotFound(ma.clone()))?;
            let mc1 = model
                .chain(mb)
                .ok_or_else(|| DockQError::ChainNotFound(mb.clone()))?;
            let small_molecule = nc0.is_het.is_some() || nc1.is_het.is_some();
            if let Some(mut info) = run_on_chains(
                (mc0, mc1),
                (nc0, nc1),
                no_align,
                capri_peptide,
                small_molecule,
            )? {
                info.chain1 = ma.clone();
                info.chain2 = mb.clone();
                result_mapping.insert(format!("{}{}", na, nb), info);
            }
        }
    }
    let total: f64 = result_mapping.values().map(|r| r.dockq).sum();
    Ok((result_mapping, total))
}

/// Parse `--mapping`. Returns (initial native→model fixed map, explicit model chain list,
/// explicit native chain list). Ports `format_mapping`.
pub(crate) fn format_mapping(
    mapping_str: &Option<String>,
) -> Result<(IndexMap<String, String>, Option<Vec<String>>, Option<Vec<String>>)> {
    let mut mapping: IndexMap<String, String> = IndexMap::new();
    let s = match mapping_str {
        Some(s) if !s.is_empty() => s,
        _ => return Ok((mapping, None, None)),
    };

    let (model_mapping, native_mapping) = match s.split_once(':') {
        Some((m, n)) => (m, n),
        None => {
            return Err(DockQError::Other(
                "When using --mapping, native chains must be set (e.g. ABC:ABC or :ABC)".into(),
            ))
        }
    };
    if native_mapping.is_empty() {
        return Err(DockQError::Other(
            "When using --mapping, native chains must be set (e.g. ABC:ABC or :ABC)".into(),
        ));
    }

    let mut model_chains = None;
    let mut native_chains = None;

    if model_mapping.is_empty() || model_mapping == "*" {
        // ":ABC" / "*:ABC" — only these natives, permute model chains.
        native_chains = Some(native_mapping.chars().map(|c| c.to_string()).collect());
    } else if model_mapping.chars().count() == native_mapping.chars().count() {
        let mm: Vec<char> = model_mapping.chars().collect();
        let nm: Vec<char> = native_mapping.chars().collect();
        for k in 0..nm.len() {
            if nm[k] != '*' && mm[k] != '*' {
                mapping.insert(nm[k].to_string(), mm[k].to_string());
            }
        }
        if *mm.last().unwrap() != '*' && *nm.last().unwrap() != '*' {
            model_chains = Some(mm.iter().map(|c| c.to_string()).collect());
            native_chains = Some(nm.iter().map(|c| c.to_string()).collect());
        }
    }
    Ok((mapping, model_chains, native_chains))
}

/// Cluster homologous model chains under each native chain (ports `group_chains`).
/// Grouping is always by sequence alignment (even under `--no_align`), counting only
/// substitution columns ('.'); a native chain with no homolog is a hard error.
fn group_chains(
    model: &Structure,
    native: &Structure,
    model_to_combo: &[String],
    native_to_combo: &[String],
    allowed_mismatches: usize,
) -> Result<(IndexMap<String, Vec<String>>, bool)> {
    let mut reverse_map = false;
    // If fewer query (model) chains than ref (native), swap (e.g. partial homomer model).
    let (qstruct, rstruct, qchains, rchains): (&Structure, &Structure, &[String], &[String]) =
        if model_to_combo.len() < native_to_combo.len() {
            reverse_map = true;
            (native, model, native_to_combo, model_to_combo)
        } else {
            (model, native, model_to_combo, native_to_combo)
        };

    let mut chain_clusters: IndexMap<String, Vec<String>> = IndexMap::new();
    for rc in rchains {
        chain_clusters.insert(rc.clone(), Vec::new());
    }

    // itertools.product(query, ref): query outer, ref inner.
    for qc_id in qchains {
        for rc_id in rchains {
            let qc = qstruct
                .chain(qc_id)
                .ok_or_else(|| DockQError::ChainNotFound(qc_id.clone()))?;
            let rc = rstruct
                .chain(rc_id)
                .ok_or_else(|| DockQError::ChainNotFound(rc_id.clone()))?;

            if qc.is_het.is_none() && rc.is_het.is_none() {
                let aln = align::align_chains(qc, rc, false);
                let n_mismatches = aln.matches.chars().filter(|&c| c == '.').count();
                if n_mismatches <= allowed_mismatches {
                    chain_clusters.get_mut(rc_id).unwrap().push(qc_id.clone());
                }
            } else if qc.is_het.is_some() && rc.is_het == qc.is_het {
                chain_clusters.get_mut(rc_id).unwrap().push(qc_id.clone());
            }
        }
    }

    let without: Vec<String> = chain_clusters
        .iter()
        .filter(|(_, v)| v.is_empty())
        .map(|(k, _)| k.clone())
        .collect();
    if !without.is_empty() {
        return Err(DockQError::NoChainMatch(without));
    }

    Ok((chain_clusters, reverse_map))
}

/// Cartesian product of pools with no element reused across positions, preserving the
/// reference enumeration order (ports `product_without_dupl`).
fn product_without_dupl(pools: &[Vec<String>]) -> Vec<Vec<String>> {
    let mut result: Vec<Vec<String>> = vec![Vec::new()];
    for pool in pools {
        let mut next: Vec<Vec<String>> = Vec::new();
        for x in &result {
            for y in pool {
                if !x.contains(y) {
                    let mut nx = x.clone();
                    nx.push(y.clone());
                    next.push(nx);
                }
            }
        }
        result = next;
    }
    result
}

/// All native→model chain maps to evaluate (ports `get_all_chain_maps`).
fn get_all_chain_maps(
    chain_clusters: &IndexMap<String, Vec<String>>,
    initial_mapping: &IndexMap<String, String>,
    reverse_map: bool,
    model_to_combo: &[String],
    native_to_combo: &[String],
) -> Vec<IndexMap<String, String>> {
    let pools: Vec<Vec<String>> = chain_clusters
        .values()
        .filter(|c| !c.is_empty())
        .cloned()
        .collect();
    let all_mappings = product_without_dupl(&pools);

    let mut maps = Vec::with_capacity(all_mappings.len());
    for mapping in &all_mappings {
        let mut chain_map = initial_mapping.clone();
        if reverse_map {
            for (i, model_chain) in model_to_combo.iter().enumerate() {
                chain_map.insert(mapping[i].clone(), model_chain.clone());
            }
        } else {
            for (i, native_chain) in native_to_combo.iter().enumerate() {
                chain_map.insert(native_chain.clone(), mapping[i].clone());
            }
        }
        maps.push(chain_map);
    }
    maps
}

/// `models:natives` string in chain-map order (ports `format_mapping_string`).
fn format_mapping_string(chain_map: &IndexMap<String, String>) -> String {
    let mut c1 = String::new();
    let mut c2 = String::new();
    for (native, model) in chain_map {
        c1.push_str(model);
        c2.push_str(native);
    }
    format!("{}:{}", c1, c2)
}

type InterfaceKey = (String, String, String, String); // (model0, model1, native0, native1)

/// Score every native interface for a single chain map, reading from the precomputed
/// cache (ports `run_on_all_native_interfaces`). Returns (per-interface results, total).
fn assemble_mapping(
    chain_map: &IndexMap<String, String>,
    cache: &HashMap<InterfaceKey, Option<InterfaceResult>>,
) -> (IndexMap<String, InterfaceResult>, f64) {
    let mut result_mapping: IndexMap<String, InterfaceResult> = IndexMap::new();
    let native_ids: Vec<&String> = chain_map.keys().collect();

    for a in 0..native_ids.len() {
        for b in (a + 1)..native_ids.len() {
            let (na, nb) = (native_ids[a], native_ids[b]);
            let (ma, mb) = (&chain_map[na], &chain_map[nb]);
            if ma == mb {
                continue; // both native chains map to the same model chain
            }
            let key = (ma.clone(), mb.clone(), na.clone(), nb.clone());
            if let Some(Some(info)) = cache.get(&key) {
                let mut info = info.clone();
                info.chain1 = ma.clone();
                info.chain2 = mb.clone();
                result_mapping.insert(format!("{}{}", na, nb), info);
            }
        }
    }
    let total: f64 = result_mapping.values().map(|r| r.dockq).sum();
    (result_mapping, total)
}

/// Score a model structure against a native structure (ports `main`'s search).
pub fn score_structures(
    model: &Structure,
    native: &Structure,
    opts: &RunOptions,
) -> Result<RunResult> {
    if opts.small_molecule {
        return Err(DockQError::SmallMoleculeUnsupported);
    }

    let (initial_mapping, model_chains_opt, native_chains_opt) = format_mapping(&opts.mapping)?;
    let model_chains: Vec<String> = model_chains_opt.unwrap_or_else(|| model.chain_ids());
    let native_chains: Vec<String> = native_chains_opt.unwrap_or_else(|| native.chain_ids());

    if model_chains.len() < 2 || native_chains.len() < 2 {
        return Err(DockQError::Other(
            "Need at least two chains in the two inputs".into(),
        ));
    }

    let fixed_models: IndexSet<&String> = initial_mapping.values().collect();
    let fixed_natives: IndexSet<&String> = initial_mapping.keys().collect();
    let model_to_combo: Vec<String> = model_chains
        .iter()
        .filter(|c| !fixed_models.contains(c))
        .cloned()
        .collect();
    let native_to_combo: Vec<String> = native_chains
        .iter()
        .filter(|c| !fixed_natives.contains(c))
        .cloned()
        .collect();

    let (chain_clusters, reverse_map) = group_chains(
        model,
        native,
        &model_to_combo,
        &native_to_combo,
        opts.allowed_mismatches,
    )?;
    let chain_maps = get_all_chain_maps(
        &chain_clusters,
        &initial_mapping,
        reverse_map,
        &model_to_combo,
        &native_to_combo,
    );
    if chain_maps.is_empty() {
        return Err(DockQError::Other(
            "No valid chain mapping could be enumerated".into(),
        ));
    }

    // Collect every distinct interface query across all maps, compute each once in parallel.
    let mut query_set: IndexSet<InterfaceKey> = IndexSet::new();
    for cm in &chain_maps {
        let nids: Vec<&String> = cm.keys().collect();
        for a in 0..nids.len() {
            for b in (a + 1)..nids.len() {
                let (na, nb) = (nids[a], nids[b]);
                let (ma, mb) = (&cm[na], &cm[nb]);
                if ma == mb {
                    continue;
                }
                query_set.insert((ma.clone(), mb.clone(), na.clone(), nb.clone()));
            }
        }
    }
    let queries: Vec<InterfaceKey> = query_set.into_iter().collect();

    let computed: Vec<(InterfaceKey, Option<InterfaceResult>)> = queries
        .par_iter()
        .map(|key| {
            let (ma, mb, na, nb) = key;
            let mc0 = model
                .chain(ma)
                .ok_or_else(|| DockQError::ChainNotFound(ma.clone()))?;
            let mc1 = model
                .chain(mb)
                .ok_or_else(|| DockQError::ChainNotFound(mb.clone()))?;
            let nc0 = native
                .chain(na)
                .ok_or_else(|| DockQError::ChainNotFound(na.clone()))?;
            let nc1 = native
                .chain(nb)
                .ok_or_else(|| DockQError::ChainNotFound(nb.clone()))?;
            let small_molecule = nc0.is_het.is_some() || nc1.is_het.is_some();
            let info = run_on_chains(
                (mc0, mc1),
                (nc0, nc1),
                opts.no_align,
                opts.capri_peptide,
                small_molecule,
            )?;
            Ok::<_, DockQError>((key.clone(), info))
        })
        .collect::<Result<Vec<_>>>()?;
    let cache: HashMap<InterfaceKey, Option<InterfaceResult>> = computed.into_iter().collect();

    // Assemble each map and pick the best by total DockQ (first wins on ties).
    let mut best_idx: Option<usize> = None;
    let mut best_total = -1.0f64;
    let mut best_result: Option<IndexMap<String, InterfaceResult>> = None;
    for (i, cm) in chain_maps.iter().enumerate() {
        let (rm, total) = assemble_mapping(cm, &cache);
        if total > best_total {
            best_total = total;
            best_idx = Some(i);
            best_result = Some(rm);
        }
    }

    let best_idx = best_idx.expect("at least one chain map scored");
    let best_result = best_result.unwrap();
    if best_result.is_empty() {
        return Err(DockQError::Other(
            "Could not find interfaces in the native. Check the inputs or select chains with --mapping.".into(),
        ));
    }
    let best_mapping = chain_maps[best_idx].clone();
    let n_interfaces = best_result.len();

    Ok(RunResult {
        model: model.id.clone(),
        native: native.id.clone(),
        best_mapping_str: format_mapping_string(&best_mapping),
        best_dockq: best_total,
        global_dockq: best_total / n_interfaces as f64,
        best_result,
        best_mapping,
    })
}
