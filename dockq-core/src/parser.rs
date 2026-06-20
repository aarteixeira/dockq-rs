//! PDB + mmCIF parsing. CONTRACT for task #3 (parser agent).
//!
//! Replicate the *observable output* of DockQ's `PDBParser` / `MMCIFParser`:
//!   - skip hydrogens (element == "H");
//!   - skip HETATM records unless `parse_hetatms`;
//!   - deduplicate altloc atoms to one-per-name (Biopython `get_atoms` representative);
//!   - build the one-letter `sequence` via seq1 (custom_map MSE->M, CME->C); for the core
//!     path het residues are skipped so the sequence is pure polymer;
//!   - mmCIF uses auth_asym_id / auth_seq_id by default;
//!   - select model `model_number` (0-based index into the file's models).
//!
//! NO SILENT FALLBACK: detect format explicitly by content (mmCIF `_atom_site` loop vs
//! PDB ATOM/HETATM records); on parse failure return `DockQError`, never warn-and-continue.
//! `pdbtbx` may be used as the loader, but the output MUST match the Python oracle dumps.

use crate::error::Result;
use crate::model::Structure;

/// Load a structure from a (optionally gzipped) PDB or mmCIF file.
///
/// * `chains` — if non-empty, restrict to these chain ids (matches DockQ's `chains=` arg).
/// * `parse_hetatms` — include HETATM records (false for the protein/NA core).
/// * `model_number` — 0-based model index to select.
pub fn load_structure(
    path: &str,
    chains: &[String],
    parse_hetatms: bool,
    model_number: usize,
) -> Result<Structure> {
    let _ = (path, chains, parse_hetatms, model_number);
    todo!("task #3: implement PDB/mmCIF parsing (Biopython-exact)")
}
