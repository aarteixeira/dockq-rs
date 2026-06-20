//! Needleman-Wunsch global alignment. CONTRACT for task #4 (alignment agent).
//!
//! Replicate Biopython's `Align.PairwiseAligner` with:
//!   match = 5, mismatch = 0, open_gap_score = -4, extend_gap_score = -0.5,
//! global mode, returning the **first** optimal alignment (`aligner.align(...)[0]`),
//! then format it exactly like DockQ's `format_alignment`:
//!   - `seq_a`: aligned model sequence (gaps as '-')
//!   - `matches`: '|' identical, ' ' if either side is '-', '.' otherwise
//!   - `seq_b`: aligned native sequence (gaps as '-')
//!
//! The first-optimal tie-break decides which residues get compared downstream, so it
//! must match Biopython. Validate against Biopython directly (it is in .venv-baseline)
//! over the example chains AND a randomized sequence-pair fuzzer.
//!
//! `use_numbering` (for --no_align): build pseudo-sequences from residue numbers
//! (each resseq mapped to chr(resn + min_resn), min_resn = max(45, -min(all_resns)))
//! and align those instead of the amino-acid sequences.

use crate::model::{Alignment, Chain};

/// Align a model chain against a native chain, returning the formatted alignment.
pub fn align_chains(model: &Chain, native: &Chain, use_numbering: bool) -> Alignment {
    let _ = (model, native, use_numbering);
    todo!("task #4: implement Biopython-compatible global alignment")
}
