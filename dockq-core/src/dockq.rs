//! Per-interface DockQ calculation (task #6, integration owner). Ports `calc_DockQ`:
//! fnat/fnonnat, iRMSD (interface backbone superposition), LRMSD (receptor-superposed
//! ligand RMSD), clashes, F1, and the DockQ formula. Small-molecule path is deferred and
//! must hard-error (`DockQError::SmallMoleculeUnsupported`), never silently fall back.

/// `DockQ = (fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(LRMSD/8.5)^2)) / 3`.
#[inline]
pub fn dockq_formula(fnat: f64, irms: f64, lrms: f64) -> f64 {
    (fnat + 1.0 / (1.0 + (irms / 1.5).powi(2)) + 1.0 / (1.0 + (lrms / 8.5).powi(2))) / 3.0
}

/// `F1 = 2*tp / (tp + fp + p)`.
#[inline]
pub fn f1(tp: f64, fp: f64, p: f64) -> f64 {
    2.0 * tp / (tp + fp + p)
}
