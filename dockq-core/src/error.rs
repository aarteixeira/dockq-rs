//! Typed errors. Policy: **no silent failures, no silent fallbacks.** Every failure
//! mode surfaces as a distinct, descriptive error rather than a warning-and-continue
//! or a quiet default.

use thiserror::Error;

#[derive(Error, Debug)]
pub enum DockQError {
    #[error("I/O error reading {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to parse structure {path}: {msg}")]
    Parse { path: String, msg: String },

    #[error("unrecognized structure format for {path}: not valid PDB or mmCIF ({detail})")]
    UnknownFormat { path: String, detail: String },

    #[error("chain '{0}' not found in the structure")]
    ChainNotFound(String),

    #[error("model index {requested} out of range (structure has {available} model(s))")]
    ModelOutOfRange { requested: usize, available: usize },

    #[error("alignment error: {0}")]
    Alignment(String),

    #[error("geometry error: {0}")]
    Geometry(String),

    #[error("native and model interfaces have incompatible sizes ({model:?} != {native:?})")]
    IncompatibleSizes {
        model: (usize, usize),
        native: (usize, usize),
    },

    #[error("no identical corresponding native chain found for: {0:?}")]
    NoChainMatch(Vec<String>),

    #[error(
        "small-molecule scoring (--small_molecule) is not implemented in this build. \
         The protein/nucleic-acid core deliberately does not silently fall back; \
         use the reference DockQ for ligand poses."
    )]
    SmallMoleculeUnsupported,

    #[error("{0}")]
    Other(String),
}

pub type Result<T> = std::result::Result<T, DockQError>;
