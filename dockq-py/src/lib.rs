//! PyO3 wrapper `dockq_rs` (task #7). Will expose drop-in `load_PDB` /
//! `run_on_all_native_interfaces`, single scoring, and the batch API.

use pyo3::prelude::*;

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn dockq_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    Ok(())
}
