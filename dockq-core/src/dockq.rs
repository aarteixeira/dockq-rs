//! Per-interface DockQ calculation. Ports `calc_DockQ` and its helpers
//! (`get_aligned_residues`, `get_residue_distances`, `get_interacting_pairs`,
//! `subset_atoms`) from the reference, using the validated geometry kernels.
//!
//! Small-molecule scoring (`calc_sym_corrected_lrmsd`) is deferred and must hard-error
//! at the call sites — never silently fall back.

use crate::constants::*;
use crate::error::{DockQError, Result};
use crate::geometry;
use crate::model::{Alignment, Chain, DistMatrix, InterfaceResult, Residue};

/// `DockQ = (fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(LRMSD/8.5)^2)) / 3`.
/// Matches the reference `dockq_formula` (which multiplies rather than `powi`, identical
/// in f64).
#[inline]
pub fn dockq_formula(fnat: f64, irms: f64, lrms: f64) -> f64 {
    (fnat + 1.0 / (1.0 + (irms / 1.5) * (irms / 1.5)) + 1.0 / (1.0 + (lrms / 8.5) * (lrms / 8.5)))
        / 3.0
}

/// `F1 = 2*tp / (tp + fp + p)`.
#[inline]
pub fn f1(tp: f64, fp: f64, p: f64) -> f64 {
    2.0 * tp / (tp + fp + p)
}

/// Aligned residue pairs from a formatted alignment (ports `get_aligned_residues`).
/// If the two aligned strings are identical, returns all residues of both chains
/// untouched (the reference's fast path). Otherwise walks the columns, advancing each
/// chain's residue cursor on non-gap, collecting pairs at '|' (exact-match) columns.
/// Iterates by `char` (not byte) so the `--no_align` numbering pseudo-sequences, which
/// may contain non-ASCII code points, align correctly.
fn get_aligned_residues<'a>(
    chain_a: &'a Chain,
    chain_b: &'a Chain,
    aln: &Alignment,
) -> (Vec<&'a Residue>, Vec<&'a Residue>) {
    if aln.seq_a == aln.seq_b {
        return (
            chain_a.residues.iter().collect(),
            chain_b.residues.iter().collect(),
        );
    }

    let sa: Vec<char> = aln.seq_a.chars().collect();
    let mt: Vec<char> = aln.matches.chars().collect();
    let sb: Vec<char> = aln.seq_b.chars().collect();

    let mut ra: Vec<&Residue> = Vec::new();
    let mut rb: Vec<&Residue> = Vec::new();
    let mut ia = 0usize;
    let mut ib = 0usize;
    let mut cur_a: Option<&Residue> = None;
    let mut cur_b: Option<&Residue> = None;

    for k in 0..sa.len() {
        if sa[k] != '-' {
            cur_a = Some(&chain_a.residues[ia]);
            ia += 1;
        }
        if sb[k] != '-' {
            cur_b = Some(&chain_b.residues[ib]);
            ib += 1;
        }
        if mt[k] == '|' {
            // At a '|' column both sides are non-gap, so both cursors are set.
            ra.push(cur_a.expect("aligned '|' column with no model residue"));
            rb.push(cur_b.expect("aligned '|' column with no native residue"));
        }
    }
    (ra, rb)
}

/// Flatten a residue group into (coords, atoms_per_res) for `residue_distances`.
/// `Residue.atoms` is already deduplicated to one-per-name, so `atoms.len()` matches the
/// reference `list_atoms_per_residue` count — keeping coord count and per-residue count
/// consistent.
fn group_coords_and_counts(group: &[&Residue]) -> (Vec<[f32; 3]>, Vec<usize>) {
    let mut coords: Vec<[f32; 3]> = Vec::new();
    let mut counts: Vec<usize> = Vec::with_capacity(group.len());
    for res in group {
        counts.push(res.atoms.len());
        for a in &res.atoms {
            coords.push(a.coord);
        }
    }
    (coords, counts)
}

/// CB coordinate, or CA if no CB (the `all_atom=False` capri-peptide path).
fn cb_or_ca(res: &Residue) -> Result<[f32; 3]> {
    if let Some(a) = res.atom_by_name("CB") {
        Ok(a.coord)
    } else if let Some(a) = res.atom_by_name("CA") {
        Ok(a.coord)
    } else {
        Err(DockQError::Geometry(format!(
            "residue {} {} has neither CB nor CA (needed for capri_peptide interface)",
            res.resname, res.resseq
        )))
    }
}

/// Residue–residue min squared-distance matrix for a residue group pair
/// (ports `get_residue_distances`). `all_atom=false` uses one CB/CA point per residue.
fn residue_group_distances(
    g1: &[&Residue],
    g2: &[&Residue],
    all_atom: bool,
) -> Result<DistMatrix> {
    if all_atom {
        let (c1, n1) = group_coords_and_counts(g1);
        let (c2, n2) = group_coords_and_counts(g2);
        Ok(geometry::residue_distances(&c1, &c2, &n1, &n2))
    } else {
        let c1: Vec<[f32; 3]> = g1.iter().map(|r| cb_or_ca(r)).collect::<Result<_>>()?;
        let c2: Vec<[f32; 3]> = g2.iter().map(|r| cb_or_ca(r)).collect::<Result<_>>()?;
        let n1 = vec![1usize; c1.len()];
        let n2 = vec![1usize; c2.len()];
        Ok(geometry::residue_distances(&c1, &c2, &n1, &n2))
    }
}

/// Unique interface residue indices on each side (ports `get_interacting_pairs` +
/// the `set(...)` dedup done by `subset_atoms`). Returns (sorted unique rows, sorted
/// unique cols) where the squared distance is `< threshold_sq`. Ordering is irrelevant
/// downstream (superposition RMSD is permutation-invariant given matched correspondence).
fn interacting_pairs(d: &DistMatrix, threshold_sq: f32) -> (Vec<usize>, Vec<usize>) {
    use std::collections::BTreeSet;
    let mut rows: BTreeSet<usize> = BTreeSet::new();
    let mut cols: BTreeSet<usize> = BTreeSet::new();
    for i in 0..d.rows {
        let base = i * d.cols;
        for j in 0..d.cols {
            if d.data[base + j] < threshold_sq {
                rows.insert(i);
                cols.insert(j);
            }
        }
    }
    (rows.into_iter().collect(), cols.into_iter().collect())
}

/// Collect paired backbone-atom coordinates (ports `subset_atoms`). For each residue
/// index (in `subset`, or all if `None`), for each atom name in `atom_types` order,
/// append (model, ref) coords only when BOTH residues have that atom. Returns
/// (model_atoms, ref_atoms) in correspondence.
fn subset_atoms(
    mod_res: &[&Residue],
    ref_res: &[&Residue],
    atom_types: &[&str],
    subset: Option<&[usize]>,
) -> (Vec<[f32; 3]>, Vec<[f32; 3]>) {
    let mut mod_atoms: Vec<[f32; 3]> = Vec::new();
    let mut ref_atoms: Vec<[f32; 3]> = Vec::new();

    let all: Vec<usize>;
    let indices: &[usize] = match subset {
        Some(s) => s,
        None => {
            all = (0..mod_res.len()).collect();
            &all
        }
    };

    for &i in indices {
        let mr = mod_res[i];
        let rr = ref_res[i];
        for &at in atom_types {
            if let (Some(ma), Some(ra)) = (mr.atom_by_name(at), rr.atom_by_name(at)) {
                mod_atoms.push(ma.coord);
                ref_atoms.push(ra.coord);
            }
        }
    }
    (mod_atoms, ref_atoms)
}

/// Superimpose `coords` onto `reference` and return the resulting RMSD
/// (mirrors `SVDSuperimposer.set(reference, coords); run(); get_rms()`).
fn superimpose_get_rms(reference: &[[f32; 3]], coords: &[[f32; 3]]) -> f32 {
    let (rot, tran) = geometry::kabsch(reference, coords);
    let transformed = geometry::apply_transform(coords, &rot, &tran);
    geometry::rmsd(&transformed, reference)
}

/// Per-interface DockQ (ports `calc_DockQ`). Returns:
///   - `Ok(None)` if the native has no contacts between the two chain groups,
///   - `Ok(Some(result))` otherwise (with `chain1`/`chain2` left blank for the caller),
///   - `Err` on a real failure (e.g. incompatible aligned sizes).
///
/// `sample_chains` / `ref_chains` are (chain group 0, chain group 1); `alignments` are
/// the model→native alignments for groups 0 and 1.
pub fn calc_dockq(
    sample_chains: (&Chain, &Chain),
    ref_chains: (&Chain, &Chain),
    alignments: (&Alignment, &Alignment),
    capri_peptide: bool,
) -> Result<Option<InterfaceResult>> {
    let fnat_threshold = if capri_peptide {
        FNAT_THRESHOLD_PEPTIDE
    } else {
        FNAT_THRESHOLD
    } as f32;
    let interface_threshold = if capri_peptide {
        INTERFACE_THRESHOLD_PEPTIDE
    } else {
        INTERFACE_THRESHOLD
    } as f32;
    let fnat_thr_sq = fnat_threshold * fnat_threshold;
    let interface_thr_sq = interface_threshold * interface_threshold;
    let clash_thr_sq = (CLASH_THRESHOLD as f32) * (CLASH_THRESHOLD as f32);

    // nat_total on the untouched native chains (full residue groups).
    let ref_all_0: Vec<&Residue> = ref_chains.0.residues.iter().collect();
    let ref_all_1: Vec<&Residue> = ref_chains.1.residues.iter().collect();
    let mut ref_res_distances = residue_group_distances(&ref_all_0, &ref_all_1, true)?;
    let nat_total = ref_res_distances.count_below(fnat_thr_sq);

    if nat_total == 0 {
        return Ok(None);
    }

    // Aligned residue subsets.
    let (aligned_sample_1, aligned_ref_1) =
        get_aligned_residues(sample_chains.0, ref_chains.0, alignments.0);
    let (aligned_sample_2, aligned_ref_2) =
        get_aligned_residues(sample_chains.1, ref_chains.1, alignments.1);

    let sample_res_distances =
        residue_group_distances(&aligned_sample_1, &aligned_sample_2, true)?;

    // If shapes differ, recompute the native matrix on the aligned native residues.
    if (ref_res_distances.rows, ref_res_distances.cols)
        != (sample_res_distances.rows, sample_res_distances.cols)
    {
        ref_res_distances = residue_group_distances(&aligned_ref_1, &aligned_ref_2, true)?;
    }

    if (sample_res_distances.rows, sample_res_distances.cols)
        != (ref_res_distances.rows, ref_res_distances.cols)
    {
        return Err(DockQError::IncompatibleSizes {
            model: (sample_res_distances.rows, sample_res_distances.cols),
            native: (ref_res_distances.rows, ref_res_distances.cols),
        });
    }

    let (nat_correct, nonnat_count, _n_native, model_total) =
        geometry::fnat_stats(&sample_res_distances, &ref_res_distances, fnat_thr_sq);

    let fnat = if nat_total != 0 {
        nat_correct as f64 / nat_total as f64
    } else {
        0.0
    };
    let fnonnat = if model_total != 0 {
        nonnat_count as f64 / model_total as f64
    } else {
        0.0
    };

    // Interface defined on the reference. For capri_peptide, recompute on CB/CA.
    if capri_peptide {
        ref_res_distances = residue_group_distances(&aligned_ref_1, &aligned_ref_2, false)?;
    }
    let (iface_rows, iface_cols) = interacting_pairs(&ref_res_distances, interface_thr_sq);

    let (sample_iface_1, ref_iface_1) = subset_atoms(
        &aligned_sample_1,
        &aligned_ref_1,
        &BACKBONE_ATOMS,
        Some(&iface_rows),
    );
    let (sample_iface_2, ref_iface_2) = subset_atoms(
        &aligned_sample_2,
        &aligned_ref_2,
        &BACKBONE_ATOMS,
        Some(&iface_cols),
    );

    let mut sample_interface_atoms = sample_iface_1;
    sample_interface_atoms.extend(sample_iface_2);
    let mut ref_interface_atoms = ref_iface_1;
    ref_interface_atoms.extend(ref_iface_2);

    // iRMSD: set(reference=sample_interface, coords=ref_interface); get_rms().
    let irms = superimpose_get_rms(&sample_interface_atoms, &ref_interface_atoms) as f64;

    // Receptor = larger native chain group (by residue count); ligand = the other.
    let ref_group1_size = ref_chains.0.len();
    let ref_group2_size = ref_chains.1.len();
    let group1_is_receptor = ref_group1_size > ref_group2_size;

    // (ref_residues, sample_residues) for receptor and ligand groups.
    let (receptor_ref, receptor_sample) = if group1_is_receptor {
        (&aligned_ref_1, &aligned_sample_1)
    } else {
        (&aligned_ref_2, &aligned_sample_2)
    };
    let (ligand_ref, ligand_sample) = if group1_is_receptor {
        (&aligned_ref_2, &aligned_sample_2)
    } else {
        (&aligned_ref_1, &aligned_sample_1)
    };
    let (class1, class2) = if group1_is_receptor {
        ("receptor", "ligand")
    } else {
        ("ligand", "receptor")
    };

    let (receptor_native, receptor_sample_atoms) =
        subset_atoms(receptor_ref, receptor_sample, &BACKBONE_ATOMS, None);
    let (ligand_native, ligand_sample_atoms) =
        subset_atoms(ligand_ref, ligand_sample, &BACKBONE_ATOMS, None);

    // LRMSD: superimpose on receptor (set(reference=receptor_native, coords=receptor_sample)),
    // apply to the ligand sample atoms, RMSD vs ligand native (no re-superposition).
    let (rot, tran) = geometry::kabsch(&receptor_native, &receptor_sample_atoms);
    let rotated_ligand = geometry::apply_transform(&ligand_sample_atoms, &rot, &tran);
    let lrms = geometry::rmsd(&ligand_native, &rotated_ligand) as f64;

    let clashes = sample_res_distances.count_below(clash_thr_sq);
    let dockq = dockq_formula(fnat, irms, lrms);
    let f1_score = f1(nat_correct as f64, nonnat_count as f64, nat_total as f64);

    Ok(Some(InterfaceResult {
        dockq,
        f1: f1_score,
        irmsd: irms,
        lrmsd: lrms,
        fnat,
        nat_correct,
        nat_total,
        fnonnat,
        nonnat_count,
        model_total,
        clashes,
        len1: ref_group1_size,
        len2: ref_group2_size,
        class1: class1.to_string(),
        class2: class2.to_string(),
        chain1: String::new(),
        chain2: String::new(),
    }))
}
