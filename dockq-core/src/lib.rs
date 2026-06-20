//! dockq-core — a from-scratch Rust reimplementation of the DockQ scoring internals
//! (protein / nucleic-acid core), built for maximum performance with a Python wrapper.
//!
//! Correctness oracle: the reference Python/Cython DockQ (v2.1.3). Every numeric kernel
//! and parser is validated against it via differential testing. No silent fallbacks.

pub mod error;
pub mod model;

pub mod align;
pub mod batch;
pub mod dockq;
pub mod geometry;
pub mod mapping;
pub mod parser;

pub use error::{DockQError, Result};
pub use model::{Alignment, Atom, Chain, DistMatrix, InterfaceResult, Residue, Structure};

pub use batch::{score_one_vs_many, score_pair, score_pairs, scan_structures, BatchOutcome};
pub use dockq::{calc_dockq, dockq_formula, f1};
pub use mapping::{run_on_native_interfaces, score_structures, RunOptions, RunResult};
pub use parser::load_structure;

/// DockQ score thresholds and atom sets (ported from `constants.py`).
pub mod constants {
    pub const FNAT_THRESHOLD: f64 = 5.0;
    pub const FNAT_THRESHOLD_PEPTIDE: f64 = 4.0;
    pub const INTERFACE_THRESHOLD: f64 = 10.0;
    pub const INTERFACE_THRESHOLD_PEPTIDE: f64 = 8.0;
    pub const CLASH_THRESHOLD: f64 = 2.0;
    pub const BOND_TOLERANCE: f64 = 0.4;

    /// Backbone atom names (protein + nucleic acid), order matters for `subset_atoms`.
    pub const BACKBONE_ATOMS: [&str; 16] = [
        "CA", "C", "N", "O", "P", "OP1", "OP2", "O2'", "O3'", "O4'", "O5'", "C1'", "C2'", "C3'",
        "C4'", "C5'",
    ];
}
