//! Batch driver. Two first-class shapes, both parallel in Rust (Rayon): one native vs
//! many models, and arbitrary (model, native) pair lists — plus a directory scan
//! convenience. Each job's parse+score runs entirely in Rust; per-job errors are returned
//! explicitly in the outcome (no silent skips). Parallelism nests safely with the
//! per-run mapping search (shared Rayon pool, no oversubscription).

use rayon::prelude::*;
use std::path::Path;

use crate::error::Result;
use crate::mapping::{format_mapping, score_structures, RunOptions, RunResult};
use crate::parser::load_structure;

/// Load + score a single model against a single native (the full single-job pipeline).
/// The `--mapping` chain filters are applied at load time, exactly as the reference does.
pub fn score_pair(model_path: &str, native_path: &str, opts: &RunOptions) -> Result<RunResult> {
    let (_initial, model_chains, native_chains) = format_mapping(&opts.mapping)?;
    let model_filter = model_chains.unwrap_or_default();
    let native_filter = native_chains.unwrap_or_default();

    let mut model = load_structure(model_path, &model_filter, opts.small_molecule, 0)?;
    let mut native = load_structure(native_path, &native_filter, opts.small_molecule, 0)?;
    // Identify structures by their source path (mirrors `model.id = path`).
    model.id = model_path.to_string();
    native.id = native_path.to_string();

    score_structures(&model, &native, opts)
}

/// Outcome of one batch job — the result or the error is carried explicitly (never dropped).
#[derive(Debug)]
pub struct BatchOutcome {
    pub model: String,
    pub native: String,
    pub result: Result<RunResult>,
}

/// Score many models against ONE native, in parallel (model-ranking / quality-assessment).
pub fn score_one_vs_many(
    native_path: &str,
    model_paths: &[String],
    opts: &RunOptions,
) -> Vec<BatchOutcome> {
    model_paths
        .par_iter()
        .map(|model_path| BatchOutcome {
            model: model_path.clone(),
            native: native_path.to_string(),
            result: score_pair(model_path, native_path, opts),
        })
        .collect()
}

/// Score arbitrary (model, native) pairs, in parallel.
pub fn score_pairs(pairs: &[(String, String)], opts: &RunOptions) -> Vec<BatchOutcome> {
    pairs
        .par_iter()
        .map(|(model_path, native_path)| BatchOutcome {
            model: model_path.clone(),
            native: native_path.clone(),
            result: score_pair(model_path, native_path, opts),
        })
        .collect()
}

/// List structure files (.pdb/.cif, optionally .gz) directly under `dir`, sorted by name.
/// Convenience for "score every model in this directory against one native".
pub fn scan_structures(dir: &str) -> Result<Vec<String>> {
    let mut out = Vec::new();
    let entries = std::fs::read_dir(dir).map_err(|e| crate::error::DockQError::Io {
        path: dir.to_string(),
        source: e,
    })?;
    for entry in entries {
        let entry = entry.map_err(|e| crate::error::DockQError::Io {
            path: dir.to_string(),
            source: e,
        })?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        if is_structure_file(&path) {
            out.push(path.to_string_lossy().into_owned());
        }
    }
    out.sort();
    Ok(out)
}

fn is_structure_file(path: &Path) -> bool {
    let name = match path.file_name().and_then(|n| n.to_str()) {
        Some(n) => n.to_ascii_lowercase(),
        None => return false,
    };
    let stem = name.strip_suffix(".gz").unwrap_or(&name);
    stem.ends_with(".pdb") || stem.ends_with(".cif") || stem.ends_with(".ent") || stem.ends_with(".mmcif")
}
