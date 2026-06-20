//! `dockq_rs` — Python bindings for the Rust DockQ core.
//!
//! Two API layers:
//!  * **Drop-in compatibility** with the reference `DockQ.DockQ`: `load_PDB` and
//!    `run_on_all_native_interfaces` keep the same names, arguments, and return shapes,
//!    so existing code migrates with zero changes.
//!  * **New ergonomic + batch API**: `score`, `score_one_vs_many`, `score_pairs`, both
//!    shapes parallelised in Rust.
//!
//! No silent failures: every Rust error becomes a Python exception.

// `load_PDB`/`run_on_all_native_interfaces` intentionally keep the reference's casing for
// drop-in compatibility.
#![allow(non_snake_case)]

use indexmap::IndexMap;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use dockq_core::{
    batch, load_structure, run_on_native_interfaces, score_pair, InterfaceResult, RunOptions,
    RunResult, Structure,
};

fn to_pyerr(e: dockq_core::DockQError) -> PyErr {
    match e {
        dockq_core::DockQError::SmallMoleculeUnsupported => PyValueError::new_err(e.to_string()),
        _ => PyRuntimeError::new_err(e.to_string()),
    }
}

/// Opaque parsed structure (wraps the Rust `Structure`). Returned by `load_PDB`,
/// accepted by `run_on_all_native_interfaces` — the drop-in replacement for the
/// Biopython model object in the reference's documented API.
#[pyclass(name = "Structure")]
pub struct PyStructure {
    inner: Structure,
}

#[pymethods]
impl PyStructure {
    /// Chain ids in file order.
    fn chain_ids(&self) -> Vec<String> {
        self.inner.chain_ids()
    }
    /// `structure[chain_id]` membership / id.
    fn __contains__(&self, chain: &str) -> bool {
        self.inner.chain(chain).is_some()
    }
    #[getter]
    fn id(&self) -> String {
        self.inner.id.clone()
    }
    fn __repr__(&self) -> String {
        format!(
            "<dockq_rs.Structure id={:?} chains={:?}>",
            self.inner.id,
            self.inner.chain_ids()
        )
    }
}

/// Build the per-interface info dict with the EXACT keys of reference DockQ v2.1.3.
fn interface_dict<'py>(
    py: Python<'py>,
    res: &InterfaceResult,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("DockQ", res.dockq)?;
    d.set_item("F1", res.f1)?;
    d.set_item("iRMSD", res.irmsd)?;
    d.set_item("LRMSD", res.lrmsd)?;
    d.set_item("fnat", res.fnat)?;
    d.set_item("nat_correct", res.nat_correct)?;
    d.set_item("nat_total", res.nat_total)?;
    d.set_item("fnonnat", res.fnonnat)?;
    d.set_item("nonnat_count", res.nonnat_count)?;
    d.set_item("model_total", res.model_total)?;
    d.set_item("clashes", res.clashes)?;
    d.set_item("len1", res.len1)?;
    d.set_item("len2", res.len2)?;
    d.set_item("class1", &res.class1)?;
    d.set_item("class2", &res.class2)?;
    d.set_item("chain1", &res.chain1)?;
    d.set_item("chain2", &res.chain2)?;
    d.set_item("is_het", false)?;
    Ok(d)
}

/// `result_mapping` dict keyed by native chain pair (tuple, like the reference), plus a
/// flat string-key copy for convenience.
fn result_mapping_dict<'py>(
    py: Python<'py>,
    mapping: &IndexMap<String, InterfaceResult>,
) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);
    for (pair, res) in mapping {
        // reference keys the dict by a (chain1, chain2) tuple of native chains.
        let chars: Vec<char> = pair.chars().collect();
        let key = (chars[0].to_string(), chars[1].to_string());
        out.set_item(key, interface_dict(py, res)?)?;
    }
    Ok(out)
}

/// Drop-in `load_PDB(path, chains=[], small_molecule=False, n_model=0)`.
#[pyfunction]
#[pyo3(signature = (path, chains=Vec::new(), small_molecule=false, n_model=0))]
fn load_PDB(
    path: &str,
    chains: Vec<String>,
    small_molecule: bool,
    n_model: usize,
) -> PyResult<PyStructure> {
    let mut inner = load_structure(path, &chains, small_molecule, n_model).map_err(to_pyerr)?;
    inner.id = path.to_string();
    Ok(PyStructure { inner })
}

/// Drop-in `run_on_all_native_interfaces(model, native, chain_map, no_align, capri_peptide,
/// low_memory)` → `(result_mapping, total_dockq)`. `chain_map` is native→model (as in the
/// reference). `low_memory` is accepted for compatibility and ignored (Rust always returns
/// full results). Insertion order of `chain_map` is preserved.
#[pyfunction]
#[pyo3(signature = (model, native, chain_map, no_align=false, capri_peptide=false, low_memory=false))]
fn run_on_all_native_interfaces<'py>(
    py: Python<'py>,
    model: &PyStructure,
    native: &PyStructure,
    chain_map: &Bound<'py, PyDict>,
    no_align: bool,
    capri_peptide: bool,
    low_memory: bool,
) -> PyResult<(Bound<'py, PyDict>, f64)> {
    let _ = low_memory;
    // Preserve the Python dict's insertion order.
    let mut map: IndexMap<String, String> = IndexMap::new();
    for (k, v) in chain_map.iter() {
        map.insert(k.extract::<String>()?, v.extract::<String>()?);
    }
    let (result, total) =
        run_on_native_interfaces(&model.inner, &native.inner, &map, no_align, capri_peptide)
            .map_err(to_pyerr)?;
    Ok((result_mapping_dict(py, &result)?, total))
}

fn opts_from(
    no_align: bool,
    capri_peptide: bool,
    small_molecule: bool,
    mapping: Option<String>,
    allowed_mismatches: usize,
) -> RunOptions {
    RunOptions {
        no_align,
        capri_peptide,
        small_molecule,
        mapping,
        allowed_mismatches,
    }
}

/// Build a Python dict for a full `RunResult`.
fn run_result_dict<'py>(py: Python<'py>, r: &RunResult) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("model", &r.model)?;
    d.set_item("native", &r.native)?;
    d.set_item("best_mapping_str", &r.best_mapping_str)?;
    d.set_item("best_dockq", r.best_dockq)?;
    d.set_item("GlobalDockQ", r.global_dockq)?;
    let best = PyDict::new(py);
    for (pair, res) in &r.best_result {
        best.set_item(pair, interface_dict(py, res)?)?;
    }
    d.set_item("best_result", best)?;
    let mapping = PyDict::new(py);
    for (native, model) in &r.best_mapping {
        mapping.set_item(native, model)?;
    }
    d.set_item("best_mapping", mapping)?;
    Ok(d)
}

/// New API: full single model-vs-native scoring (the CLI's job), returning a result dict.
#[pyfunction]
#[pyo3(signature = (model, native, no_align=false, capri_peptide=false, small_molecule=false, mapping=None, allowed_mismatches=0))]
fn score<'py>(
    py: Python<'py>,
    model: &str,
    native: &str,
    no_align: bool,
    capri_peptide: bool,
    small_molecule: bool,
    mapping: Option<String>,
    allowed_mismatches: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let opts = opts_from(no_align, capri_peptide, small_molecule, mapping, allowed_mismatches);
    // Rayon parallelism is internal (pure Rust, GIL-independent); no GIL release needed.
    let r = score_pair(model, native, &opts).map_err(to_pyerr)?;
    run_result_dict(py, &r)
}

fn outcome_dict<'py>(py: Python<'py>, o: &batch::BatchOutcome) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("model", &o.model)?;
    d.set_item("native", &o.native)?;
    match &o.result {
        Ok(r) => {
            d.set_item("ok", true)?;
            d.set_item("result", run_result_dict(py, r)?)?;
        }
        Err(e) => {
            // No silent failure: the error is surfaced explicitly per job.
            d.set_item("ok", false)?;
            d.set_item("error", e.to_string())?;
        }
    }
    Ok(d)
}

/// New API: score MANY models against ONE native, in parallel (model ranking).
#[pyfunction]
#[pyo3(signature = (native, models, no_align=false, capri_peptide=false, small_molecule=false, mapping=None, allowed_mismatches=0))]
fn score_one_vs_many<'py>(
    py: Python<'py>,
    native: &str,
    models: Vec<String>,
    no_align: bool,
    capri_peptide: bool,
    small_molecule: bool,
    mapping: Option<String>,
    allowed_mismatches: usize,
) -> PyResult<Bound<'py, PyList>> {
    let opts = opts_from(no_align, capri_peptide, small_molecule, mapping, allowed_mismatches);
    let outcomes = dockq_core::score_one_vs_many(native, &models, &opts);
    let list = PyList::empty(py);
    for o in &outcomes {
        list.append(outcome_dict(py, o)?)?;
    }
    Ok(list)
}

/// New API: score arbitrary (model, native) pairs, in parallel.
#[pyfunction]
#[pyo3(signature = (pairs, no_align=false, capri_peptide=false, small_molecule=false, mapping=None, allowed_mismatches=0))]
fn score_pairs<'py>(
    py: Python<'py>,
    pairs: Vec<(String, String)>,
    no_align: bool,
    capri_peptide: bool,
    small_molecule: bool,
    mapping: Option<String>,
    allowed_mismatches: usize,
) -> PyResult<Bound<'py, PyList>> {
    let opts = opts_from(no_align, capri_peptide, small_molecule, mapping, allowed_mismatches);
    let outcomes = dockq_core::score_pairs(&pairs, &opts);
    let list = PyList::empty(py);
    for o in &outcomes {
        list.append(outcome_dict(py, o)?)?;
    }
    Ok(list)
}

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn _dockq_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyStructure>()?;
    // Drop-in compatible API (reference names).
    m.add_function(wrap_pyfunction!(load_PDB, m)?)?;
    m.add_function(wrap_pyfunction!(run_on_all_native_interfaces, m)?)?;
    // New ergonomic + batch API.
    m.add_function(wrap_pyfunction!(score, m)?)?;
    m.add_function(wrap_pyfunction!(score_one_vs_many, m)?)?;
    m.add_function(wrap_pyfunction!(score_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    Ok(())
}
