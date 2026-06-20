//! Needleman-Wunsch global alignment. CONTRACT task #4 (alignment agent).
//!
//! Replicate Biopython's `Align.PairwiseAligner` with:
//! match = 5, mismatch = 0, open_gap_score = -4, extend_gap_score = -0.5,
//! global mode, returning **first** optimal alignment (`aligner.align(...)[0]`),
//! format exactly DockQ's `format_alignment`:
//! - `seq_a`: aligned model sequence (gaps '-')
//! - `matches`: '|' identical, side '-', '.' otherwise
//! - `seq_b`: aligned native sequence (gaps as '-')
//!
//! first-optimal tie-break decides residues get compared downstream,
//! match Biopython. Validate Biopython directly in .venv-baseline)
//! over example chains AND randomized sequence-pair fuzzer.
//!
//! `use_numbering` --no_align): build pseudo-sequences residue numbers
//! resseq mapped chr(resn + min_resn), min_resn = max(45, -min(all_resns)))
//! align instead amino-acid sequences.
//!
//! ## Reverse-engineered Biopython 1.87 details (grounded in
//! `Bio/Align/_pairwisealigner.c`, function `PathGenerator_next_gotoh_global`,
//! and confirmed empirically against the live `.venv-baseline`):
//!
//! * **Affine gap convention (Gotoh):** a gap of length `L` scores
//!   `open + (L-1)*extend`. Entering a gap from the match state costs `open`
//!   (= -4); each additional gap column costs `extend` (= -0.5).
//! * **End gaps** are scored exactly like internal gaps (DockQ leaves
//!   `*_end_*` scores at default, which equals the internal gap score).
//! * **First-optimal tie-break:** Biopython's `[0]` traceback walks from the
//!   bottom-right corner and, at every cell, prefers predecessor matrices in
//!   the order **M (diagonal) > Ix (vertical, model-residue vs gap) >
//!   Iy (horizontal, gap vs native-residue)** — see the `while(1)` loop at
//!   lines ~1030-1049 of `_pairwisealigner.c` which tests `M_MATRIX`, then
//!   `Ix_MATRIX`, then `Iy_MATRIX`. Because diagonal moves are consumed first
//!   while tracing back (right-to-left), this pushes gaps as far toward the
//!   **start** of the alignment as the optimal score allows.
//!   - `VERTICAL` advances the target/model index `i` only -> model char in
//!     `seq_a`, '-' in `seq_b`.
//!   - `HORIZONTAL` advances the query/native index `j` only -> '-' in
//!     `seq_a`, native char in `seq_b`.

use crate::model::{Alignment, Chain};

const MATCH: f64 = 5.0;
const MISMATCH: f64 = 0.0;
const OPEN: f64 = -4.0;
const EXTEND: f64 = -0.5;

/// Tolerance for floating-point equality when collecting optimal predecessors.
/// Biopython compares the C `double` accumulators exactly; our DP performs the
/// same additions in the same order, so an exact `==` would in principle work,
/// but a tiny epsilon guards against benign reassociation without ever merging
/// genuinely distinct scores (the score lattice here is spaced by 0.5).
const EPS: f64 = 1e-9;

/// Predecessor bitmask for a DP cell. Mirrors Biopython's trace bits
/// (`M_MATRIX = 0x1`, `Ix_MATRIX = 0x2`, `Iy_MATRIX = 0x4`).
const FROM_M: u8 = 0x1;
const FROM_IX: u8 = 0x2;
const FROM_IY: u8 = 0x4;

/// Align model chain native chain, returning formatted alignment.
pub fn align_chains(model: &Chain, native: &Chain, use_numbering: bool) -> Alignment {
    let (seq_a_src, seq_b_src): (Vec<char>, Vec<char>) = if use_numbering {
        numbering_pseudo_sequences(model, native)
    } else {
        (
            model.sequence.chars().collect(),
            native.sequence.chars().collect(),
        )
    };
    align_seqs(&seq_a_src, &seq_b_src)
}

/// Build the residue-number pseudo-sequences for `--no_align` mode, exactly as
/// DockQ's `align_chains(..., use_numbering=True)`:
///
/// ```text
/// min_resn = max(45, -min(model_resns + native_resns))
/// seq      = ''.join(chr(resn + min_resn) for resn in numbering)
/// ```
///
/// `chr()` can exceed ASCII; Python `str` holds Unicode code points, so we use
/// `char::from_u32` to mirror that. Invalid code points (surrogates / > 0x10FFFF)
/// cannot arise for realistic resseq values but, if they did, Python's `chr`
/// would raise; here we fall back to the replacement character so the function
/// stays total. (resseq is `i64` in our model; DockQ uses `int(residue.id[1])`.)
fn numbering_pseudo_sequences(model: &Chain, native: &Chain) -> (Vec<char>, Vec<char>) {
    let model_nums: Vec<i64> = model.residues.iter().map(|r| r.resseq).collect();
    let native_nums: Vec<i64> = native.residues.iter().map(|r| r.resseq).collect();

    let min_all = model_nums
        .iter()
        .chain(native_nums.iter())
        .copied()
        .min()
        .unwrap_or(0);
    // max(45, -min(...))
    let min_resn: i64 = std::cmp::max(45, -min_all);

    let to_char = |resn: i64| -> char {
        let code = resn + min_resn;
        u32::try_from(code)
            .ok()
            .and_then(char::from_u32)
            .unwrap_or('\u{FFFD}')
    };

    (
        model_nums.iter().map(|&n| to_char(n)).collect(),
        native_nums.iter().map(|&n| to_char(n)).collect(),
    )
}

/// Score of aligning two characters (the substitution score).
#[inline]
fn sub_score(a: char, b: char) -> f64 {
    if a == b {
        MATCH
    } else {
        MISMATCH
    }
}

/// Core global affine alignment + Biopython-compatible first-path traceback.
///
/// `a` is the target (model, becomes `seq_a`); `b` is the query (native,
/// becomes `seq_b`).
fn align_seqs(a: &[char], b: &[char]) -> Alignment {
    let n = a.len();
    let m = b.len();

    // Degenerate cases: if either side is empty, the only alignment is all gaps.
    if n == 0 && m == 0 {
        return Alignment {
            seq_a: String::new(),
            matches: String::new(),
            seq_b: String::new(),
        };
    }

    const NEG_INF: f64 = f64::NEG_INFINITY;

    // Three score matrices of size (n+1) x (m+1), row-major.
    // m_mat[i][j]: best score of a..i vs b..j ending with i,j aligned (diagonal).
    // ix_mat[i][j]: best score ending with a gap in the QUERY (vertical: a[i-1] vs '-').
    // iy_mat[i][j]: best score ending with a gap in the TARGET (horizontal: '-' vs b[j-1]).
    let w = m + 1;
    let idx = |i: usize, j: usize| i * w + j;

    let mut m_mat = vec![NEG_INF; (n + 1) * w];
    let mut ix_mat = vec![NEG_INF; (n + 1) * w];
    let mut iy_mat = vec![NEG_INF; (n + 1) * w];

    // Trace bitmasks: which predecessor matrices achieve the optimum for this cell.
    let mut m_tr = vec![0u8; (n + 1) * w];
    let mut ix_tr = vec![0u8; (n + 1) * w];
    let mut iy_tr = vec![0u8; (n + 1) * w];

    // --- Initialization (corner + first row/col) ---
    m_mat[idx(0, 0)] = 0.0;

    // First column (j = 0): only vertical gaps (consume target residues a[0..i]).
    // ix[i][0] = gap of length i in query = open + (i-1)*extend.
    for i in 1..=n {
        // From M (i==1) -> open; from Ix (i>1) -> extend.
        let from_m = add(m_mat[idx(i - 1, 0)], OPEN);
        let from_ix = add(ix_mat[idx(i - 1, 0)], EXTEND);
        let (best, tr) = best_gap(from_m, from_ix);
        ix_mat[idx(i, 0)] = best;
        ix_tr[idx(i, 0)] = tr;
    }

    // First row (i = 0): only horizontal gaps (consume query residues b[0..j]).
    for j in 1..=m {
        let from_m = add(m_mat[idx(0, j - 1)], OPEN);
        let from_iy = add(iy_mat[idx(0, j - 1)], EXTEND);
        let (best, tr) = best_gap_y(from_m, from_iy);
        iy_mat[idx(0, j)] = best;
        iy_tr[idx(0, j)] = tr;
    }

    // --- Main DP fill ---
    for i in 1..=n {
        for j in 1..=m {
            // M[i][j]: align a[i-1] with b[j-1]; predecessor is best of the three
            // matrices at (i-1, j-1), plus the substitution score.
            {
                let s = sub_score(a[i - 1], b[j - 1]);
                let pm = m_mat[idx(i - 1, j - 1)];
                let px = ix_mat[idx(i - 1, j - 1)];
                let py = iy_mat[idx(i - 1, j - 1)];
                let (best, tr) = best3(pm, px, py);
                m_mat[idx(i, j)] = if best > NEG_INF { best + s } else { NEG_INF };
                m_tr[idx(i, j)] = tr;
            }

            // Ix[i][j]: gap in query (vertical). Come from (i-1, j): from M -> open,
            // from Ix -> extend. (Iy -> Ix would also be an "open"; Biopython's
            // gotoh allows M->Ix and Ix->Ix; transitions from Iy into Ix are NOT
            // permitted in the standard 3-state gotoh, matching PairwiseAligner.)
            {
                let from_m = add(m_mat[idx(i - 1, j)], OPEN);
                let from_ix = add(ix_mat[idx(i - 1, j)], EXTEND);
                let (best, tr) = best_gap(from_m, from_ix);
                ix_mat[idx(i, j)] = best;
                ix_tr[idx(i, j)] = tr;
            }

            // Iy[i][j]: gap in target (horizontal). Come from (i, j-1): from M -> open,
            // from Iy -> extend.
            {
                let from_m = add(m_mat[idx(i, j - 1)], OPEN);
                let from_iy = add(iy_mat[idx(i, j - 1)], EXTEND);
                let (best, tr) = best_gap_y(from_m, from_iy);
                iy_mat[idx(i, j)] = best;
                iy_tr[idx(i, j)] = tr;
            }
        }
    }

    // --- Determine the starting matrix at the bottom-right corner ---
    // Biopython's first path prefers M, then Ix, then Iy among matrices whose
    // corner score equals the global optimum.
    let corner_m = m_mat[idx(n, m)];
    let corner_ix = ix_mat[idx(n, m)];
    let corner_iy = iy_mat[idx(n, m)];
    let best_corner = corner_m.max(corner_ix).max(corner_iy);

    #[derive(Clone, Copy, PartialEq)]
    enum Mat {
        M,
        Ix,
        Iy,
    }

    let mut state = if (corner_m - best_corner).abs() <= EPS {
        Mat::M
    } else if (corner_ix - best_corner).abs() <= EPS {
        Mat::Ix
    } else {
        Mat::Iy
    };

    // --- Traceback (first path) ---
    // We emit columns from the bottom-right toward the top-left, then reverse.
    // At each step, the current matrix tells us the move:
    //   M  -> DIAGONAL: a[i-1] vs b[j-1]; i--, j--
    //   Ix -> VERTICAL: a[i-1] vs '-';    i--
    //   Iy -> HORIZONTAL: '-' vs b[j-1];  j--
    // The next matrix is chosen from the current cell's trace bits, preferring
    // M > Ix > Iy (Biopython's order).
    let mut i = n;
    let mut j = m;
    let mut col_a: Vec<char> = Vec::with_capacity(n + m);
    let mut col_b: Vec<char> = Vec::with_capacity(n + m);

    let pick_next = |tr: u8| -> Mat {
        if tr & FROM_M != 0 {
            Mat::M
        } else if tr & FROM_IX != 0 {
            Mat::Ix
        } else {
            Mat::Iy
        }
    };

    while i > 0 || j > 0 {
        match state {
            Mat::M => {
                // diagonal
                col_a.push(a[i - 1]);
                col_b.push(b[j - 1]);
                let tr = m_tr[idx(i, j)];
                i -= 1;
                j -= 1;
                state = pick_next(tr);
            }
            Mat::Ix => {
                // vertical: model char, gap in native
                col_a.push(a[i - 1]);
                col_b.push('-');
                let tr = ix_tr[idx(i, j)];
                i -= 1;
                state = pick_next(tr);
            }
            Mat::Iy => {
                // horizontal: gap in model, native char
                col_a.push('-');
                col_b.push(b[j - 1]);
                let tr = iy_tr[idx(i, j)];
                j -= 1;
                state = pick_next(tr);
            }
        }
    }

    col_a.reverse();
    col_b.reverse();

    let seq_a: String = col_a.iter().collect();
    let seq_b: String = col_b.iter().collect();
    let matches: String = col_a
        .iter()
        .zip(col_b.iter())
        .map(|(&x, &y)| {
            if x == y {
                '|'
            } else if x == '-' || y == '-' {
                ' '
            } else {
                '.'
            }
        })
        .collect();

    Alignment {
        seq_a,
        matches,
        seq_b,
    }
}

/// Add with -inf absorption (so -inf + finite stays -inf, never NaN).
#[inline]
fn add(x: f64, d: f64) -> f64 {
    if x == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    } else {
        x + d
    }
}

/// Best of three predecessor scores (for the M matrix), returning the score and
/// the set of optimal-predecessor bits (M / Ix / Iy), preferring none — all
/// equal-optimum predecessors are recorded; priority is applied later during
/// traceback.
#[inline]
fn best3(pm: f64, px: f64, py: f64) -> (f64, u8) {
    let best = pm.max(px).max(py);
    if best == f64::NEG_INFINITY {
        return (f64::NEG_INFINITY, 0);
    }
    let mut tr = 0u8;
    if (pm - best).abs() <= EPS && pm > f64::NEG_INFINITY {
        tr |= FROM_M;
    }
    if (px - best).abs() <= EPS && px > f64::NEG_INFINITY {
        tr |= FROM_IX;
    }
    if (py - best).abs() <= EPS && py > f64::NEG_INFINITY {
        tr |= FROM_IY;
    }
    (best, tr)
}

/// Best predecessor for the Ix (vertical-gap) matrix: from M (open) or Ix (extend).
#[inline]
fn best_gap(from_m: f64, from_ix: f64) -> (f64, u8) {
    let best = from_m.max(from_ix);
    if best == f64::NEG_INFINITY {
        return (f64::NEG_INFINITY, 0);
    }
    let mut tr = 0u8;
    if (from_m - best).abs() <= EPS && from_m > f64::NEG_INFINITY {
        tr |= FROM_M;
    }
    if (from_ix - best).abs() <= EPS && from_ix > f64::NEG_INFINITY {
        tr |= FROM_IX;
    }
    (best, tr)
}

/// Best predecessor for the Iy (horizontal-gap) matrix: from M (open) or Iy (extend).
#[inline]
fn best_gap_y(from_m: f64, from_iy: f64) -> (f64, u8) {
    let best = from_m.max(from_iy);
    if best == f64::NEG_INFINITY {
        return (f64::NEG_INFINITY, 0);
    }
    let mut tr = 0u8;
    if (from_m - best).abs() <= EPS && from_m > f64::NEG_INFINITY {
        tr |= FROM_M;
    }
    if (from_iy - best).abs() <= EPS && from_iy > f64::NEG_INFINITY {
        tr |= FROM_IY;
    }
    (best, tr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Chain, Residue};

    /// Build a polymer Chain from a one-letter sequence (residue numbers 1..=len).
    fn chain_from_seq(id: &str, seq: &str) -> Chain {
        let residues = seq
            .chars()
            .enumerate()
            .map(|(k, c)| Residue {
                het_flag: ' ',
                resseq: (k + 1) as i64,
                icode: ' ',
                resname: c.to_string(),
                resname1: c.to_string(),
                atoms: Vec::new(),
            })
            .collect();
        Chain {
            id: id.to_string(),
            residues,
            sequence: seq.to_string(),
            is_het: None,
        }
    }

    /// Build a Chain with explicit residue numbers (for use_numbering tests).
    fn chain_from_nums(id: &str, nums: &[i64]) -> Chain {
        let residues = nums
            .iter()
            .map(|&n| Residue {
                het_flag: ' ',
                resseq: n,
                icode: ' ',
                resname: "ALA".to_string(),
                resname1: "A".to_string(),
                atoms: Vec::new(),
            })
            .collect();
        // sequence is irrelevant in numbering mode; fill with A's.
        Chain {
            id: id.to_string(),
            residues,
            sequence: "A".repeat(nums.len()),
            is_het: None,
        }
    }

    fn align(a: &str, b: &str) -> (String, String, String) {
        let ca = chain_from_seq("A", a);
        let cb = chain_from_seq("B", b);
        let al = align_chains(&ca, &cb, false);
        (al.seq_a, al.matches, al.seq_b)
    }

    #[test]
    fn identical_sequences_all_match_no_gaps() {
        let (a, m, b) = align("ACDEFGHIK", "ACDEFGHIK");
        assert_eq!(a, "ACDEFGHIK");
        assert_eq!(m, "|||||||||");
        assert_eq!(b, "ACDEFGHIK");
    }

    #[test]
    fn single_substitution_gives_dot() {
        // One mismatch in the middle. Biopython: no gaps (mismatch=0 still aligns).
        let (a, m, b) = align("ACDEF", "ACXEF");
        assert_eq!(a, "ACDEF");
        assert_eq!(m, "||.||");
        assert_eq!(b, "ACXEF");
    }

    #[test]
    fn deletion_inserts_gap_in_query() {
        // model "AA" vs native "A": single deletion. Biopython [0] -> gap at LEFT
        // (in seq_b). Confirmed live: seqA="AA", seqB="-A".
        let (a, m, b) = align("AA", "A");
        assert_eq!(a, "AA");
        assert_eq!(m, " |");
        assert_eq!(b, "-A");
    }

    #[test]
    fn insertion_inserts_gap_in_target() {
        // model "A" vs native "AA": Biopython [0] -> seqA="-A", seqB="AA".
        let (a, m, b) = align("A", "AA");
        assert_eq!(a, "-A");
        assert_eq!(m, " |");
        assert_eq!(b, "AA");
    }

    #[test]
    fn tie_break_repeats_gaps_pushed_left() {
        // "AAAA" vs "AA": 3 optimal alignments tie; Biopython [0] -> "--AA".
        // Proves our gaps-to-the-left tie-break matches.
        let (a, m, b) = align("AAAA", "AA");
        assert_eq!(a, "AAAA");
        assert_eq!(m, "  ||");
        assert_eq!(b, "--AA");
    }

    #[test]
    fn tie_break_leading_deletion_choice() {
        // "AAB" vs "AB": which A is deleted? Biopython [0] -> seqB="-AB".
        let (a, m, b) = align("AAB", "AB");
        assert_eq!(a, "AAB");
        assert_eq!(m, " ||");
        assert_eq!(b, "-AB");
    }

    #[test]
    fn real_example_chain_pair_model_native_a() {
        // The model/native example A chains differ at one position (E[I/P]QRTPK).
        // Full sequences from examples/model.pdb and native.pdb (chain A).
        let model_a = "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM";
        let native_a = "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM";
        let (a, m, b) = align(model_a, native_a);
        // No gaps: equal length, one mismatch.
        assert_eq!(a, model_a);
        assert_eq!(b, native_a);
        assert_eq!(a.len(), m.len());
        assert_eq!(b.len(), m.len());
        // Exactly one '.' (the I/P difference), rest '|', no spaces.
        assert_eq!(m.matches('.').count(), 1);
        assert_eq!(m.matches(' ').count(), 0);
        assert_eq!(m.matches('|').count(), model_a.len() - 1);
    }

    #[test]
    fn use_numbering_simple() {
        // Numbering mode: model resns [1,2,3], native [1,2,3] -> identical pseudo-seq.
        let cm = chain_from_nums("A", &[1, 2, 3]);
        let cn = chain_from_nums("B", &[1, 2, 3]);
        let al = align_chains(&cm, &cn, true);
        assert_eq!(al.matches, "|||");
        assert_eq!(al.seq_a.chars().count(), 3);
        assert_eq!(al.seq_b.chars().count(), 3);
    }

    #[test]
    fn use_numbering_with_gap() {
        // model resns [1,2,3,4], native [1,2,4] -> native is missing residue 3.
        // Pseudo-seqs (min_resn=45): model = chars(46,47,48,49), native=(46,47,49).
        // Best alignment deletes the '3' position. Verify a single internal gap.
        let cm = chain_from_nums("A", &[1, 2, 3, 4]);
        let cn = chain_from_nums("B", &[1, 2, 4]);
        let al = align_chains(&cm, &cn, true);
        assert_eq!(al.seq_a.chars().count(), 4);
        assert_eq!(al.seq_b.chars().count(), 4);
        // One gap in seq_b.
        assert_eq!(al.seq_b.matches('-').count(), 1);
        assert_eq!(al.matches.matches('|').count(), 3);
    }

    #[test]
    fn empty_sequences() {
        let (a, m, b) = align("", "");
        assert_eq!(a, "");
        assert_eq!(m, "");
        assert_eq!(b, "");
    }

    /// Differential fuzzer against live Biopython output dumped to
    /// `/tmp/align_oracle.json` by `/tmp/gen_align_oracle.py` (see report).
    /// Each case provides Biopython's `[0]` formatted alignment; we require an
    /// EXACT string match of seqA / matches / seqB. Ignored by default (needs
    /// the external JSON); run with `cargo test -p dockq-core -- --ignored fuzz`.
    #[test]
    #[ignore]
    fn fuzz_against_biopython_oracle() {
        use std::collections::BTreeMap;
        let path = std::env::var("DOCKQ_ALIGN_ORACLE").unwrap_or_else(|_| {
            std::env::temp_dir()
                .join("align_oracle.json")
                .to_string_lossy()
                .into_owned()
        });
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {path}: {e}"));

        // Minimal hand-rolled extraction of the JSON string fields we need,
        // avoiding any new dependency. The generator emits a flat array under
        // "cases" with string fields a,b,seqA,matches,seqB,cat. We parse with a
        // tiny state machine over the serde_json-style escaping the generator
        // produced (only \" and \\ and unicode are possible; our sequences are
        // ASCII letters / '-' so no escapes appear in seq fields, but pseudo
        // numbering chars never enter this path).
        let cases = parse_cases(&data);
        assert!(!cases.is_empty(), "no cases parsed from oracle");

        let mut fails: Vec<String> = Vec::new();
        let mut per_cat_total: BTreeMap<String, usize> = BTreeMap::new();
        let mut per_cat_fail: BTreeMap<String, usize> = BTreeMap::new();

        for c in &cases {
            *per_cat_total.entry(c.cat.clone()).or_default() += 1;
            let a_chars: Vec<char> = c.a.chars().collect();
            let b_chars: Vec<char> = c.b.chars().collect();
            let got = super::align_seqs(&a_chars, &b_chars);
            if got.seq_a != c.seq_a || got.matches != c.matches || got.seq_b != c.seq_b {
                *per_cat_fail.entry(c.cat.clone()).or_default() += 1;
                if fails.len() < 20 {
                    fails.push(format!(
                        "cat={} a={:?} b={:?}\n  got  A={:?}\n  want A={:?}\n  got  M={:?}\n  want M={:?}\n  got  B={:?}\n  want B={:?}",
                        c.cat, c.a, c.b, got.seq_a, c.seq_a, got.matches, c.matches, got.seq_b, c.seq_b
                    ));
                }
            }
        }

        let total = cases.len();
        let failed: usize = per_cat_fail.values().sum();
        eprintln!("FUZZ: {} cases, {} failed", total, failed);
        for (cat, n) in &per_cat_total {
            let f = per_cat_fail.get(cat).copied().unwrap_or(0);
            eprintln!("  {cat}: {n} cases, {f} fail");
        }
        if failed > 0 {
            for f in &fails {
                eprintln!("--- MISMATCH ---\n{f}");
            }
            panic!("{failed}/{total} fuzzer cases mismatched Biopython");
        }
    }

    /// Differential fuzzer for `use_numbering` mode against live DockQ
    /// `align_chains(..., use_numbering=True)`, dumped to
    /// `/tmp/numbering_oracle.json` by `/tmp/gen_numbering_oracle.py`. Each case
    /// gives the model/native resseq lists and the expected formatted alignment;
    /// we rebuild chains, run `align_chains(.., true)`, and require an EXACT match.
    #[test]
    #[ignore]
    fn fuzz_numbering_against_dockq() {
        let path = std::env::var("DOCKQ_NUMBERING_ORACLE").unwrap_or_else(|_| {
            std::env::temp_dir()
                .join("numbering_oracle.json")
                .to_string_lossy()
                .into_owned()
        });
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {path}: {e}"));
        let cases = parse_numbering_cases(&data);
        assert!(!cases.is_empty(), "no numbering cases parsed");

        let mut failed = 0usize;
        let mut shown = 0usize;
        for c in &cases {
            let cm = chain_from_nums("M", &c.model_nums);
            let cn = chain_from_nums("N", &c.native_nums);
            let got = align_chains(&cm, &cn, true);
            if got.seq_a != c.seq_a || got.matches != c.matches || got.seq_b != c.seq_b {
                failed += 1;
                if shown < 10 {
                    shown += 1;
                    eprintln!(
                        "NUM MISMATCH {}\n  got  A={:?}\n  want A={:?}\n  got  M={:?}\n  want M={:?}\n  got  B={:?}\n  want B={:?}",
                        c.cat, got.seq_a, c.seq_a, got.matches, c.matches, got.seq_b, c.seq_b
                    );
                }
            }
        }
        eprintln!("NUMBERING FUZZ: {} cases, {} failed", cases.len(), failed);
        assert_eq!(failed, 0, "{failed}/{} numbering cases mismatched", cases.len());
    }

    struct NumCase {
        cat: String,
        model_nums: Vec<i64>,
        native_nums: Vec<i64>,
        seq_a: String,
        matches: String,
        seq_b: String,
    }

    fn parse_numbering_cases(data: &str) -> Vec<NumCase> {
        split_case_objects(data)
            .into_iter()
            .map(|obj| NumCase {
                cat: field(obj, "cat").unwrap_or_default(),
                model_nums: int_array(obj, "model_nums"),
                native_nums: int_array(obj, "native_nums"),
                seq_a: field(obj, "seqA").unwrap_or_default(),
                matches: field(obj, "matches").unwrap_or_default(),
                seq_b: field(obj, "seqB").unwrap_or_default(),
            })
            .collect()
    }

    /// Split the `"cases": [ {..}, {..}, ... ]` array into the slices of each
    /// top-level object. STRING-AWARE: brackets/braces inside JSON string values
    /// are ignored (the numbering pseudo-sequences contain literal `[`, `]`,
    /// `{`, `}`, `\` characters from `chr()`, which would otherwise corrupt a
    /// naive brace counter).
    fn split_case_objects(data: &str) -> Vec<&str> {
        let key = "\"cases\":";
        let start = match data.find(key) {
            Some(s) => s + key.len(),
            None => return Vec::new(),
        };
        let arr = &data[start..];
        let arr_start = match arr.find('[') {
            Some(p) => p,
            None => return Vec::new(),
        };
        let bytes = arr.as_bytes();
        let mut out = Vec::new();
        let mut i = arr_start + 1;
        let mut depth = 1usize; // inside the cases array
        let mut obj_start: Option<usize> = None;
        let mut in_str = false;
        let mut escaped = false;
        while i < bytes.len() {
            let ch = bytes[i];
            if in_str {
                if escaped {
                    escaped = false;
                } else if ch == b'\\' {
                    escaped = true;
                } else if ch == b'"' {
                    in_str = false;
                }
                i += 1;
                continue;
            }
            match ch {
                b'"' => in_str = true,
                b'{' => {
                    if depth == 1 {
                        obj_start = Some(i);
                    }
                    depth += 1;
                }
                b'}' => {
                    depth -= 1;
                    if depth == 1 {
                        if let Some(s) = obj_start.take() {
                            out.push(&arr[s..=i]);
                        }
                    }
                }
                b'[' => depth += 1,
                b']' => {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                }
                _ => {}
            }
            i += 1;
        }
        out
    }

    /// Extract `"name": [n, n, ...]` integer array from a flat JSON object slice.
    fn int_array(obj: &str, name: &str) -> Vec<i64> {
        let pat = format!("\"{name}\":");
        let p = match obj.find(&pat) {
            Some(p) => p + pat.len(),
            None => return Vec::new(),
        };
        let rest = obj[p..].trim_start();
        let bytes = rest.as_bytes();
        if bytes.first() != Some(&b'[') {
            return Vec::new();
        }
        let end = rest.find(']').unwrap_or(rest.len());
        rest[1..end]
            .split(',')
            .filter_map(|s| {
                let t = s.trim();
                if t.is_empty() {
                    None
                } else {
                    t.parse::<i64>().ok()
                }
            })
            .collect()
    }

    /// A single fuzzer case.
    struct Case {
        cat: String,
        a: String,
        b: String,
        seq_a: String,
        matches: String,
        seq_b: String,
    }

    /// Extract the `cases` array from the oracle JSON. The generator uses
    /// `json.dump` (compact), so each object looks like
    /// `{"cat": "...", "a": "...", "b": "...", "seqA": "...", "matches": "...", "seqB": "..."}`.
    /// Sequence/match strings contain only ASCII letters and '-' (no escapes),
    /// so a tolerant string-field scan suffices and avoids a JSON dependency.
    fn parse_cases(data: &str) -> Vec<Case> {
        split_case_objects(data)
            .into_iter()
            .filter_map(parse_obj)
            .collect()
    }

    fn parse_obj(obj: &str) -> Option<Case> {
        Some(Case {
            cat: field(obj, "cat")?,
            a: field(obj, "a")?,
            b: field(obj, "b")?,
            seq_a: field(obj, "seqA")?,
            matches: field(obj, "matches")?,
            seq_b: field(obj, "seqB")?,
        })
    }

    /// Extract `"name": "value"` string value from a flat JSON object slice,
    /// honoring `\"` and `\\` escapes (the only ones serde/json.dump emit for
    /// our content). Returns the unescaped value.
    fn field(obj: &str, name: &str) -> Option<String> {
        let pat = format!("\"{name}\":");
        let mut search = 0usize;
        // The keys are unique within the object, but "a" is a prefix of nothing
        // problematic because we anchor on the full quoted key + colon.
        let p = obj[search..].find(&pat)? + search;
        search = p + pat.len();
        let rest = obj[search..].trim_start();
        let rb = rest.as_bytes();
        if rb.first() != Some(&b'"') {
            return None;
        }
        let mut val = String::new();
        let mut k = 1usize; // skip opening quote
        let bytes = rest.as_bytes();
        while k < bytes.len() {
            let ch = bytes[k];
            if ch == b'\\' {
                k += 1;
                if k >= bytes.len() {
                    break;
                }
                match bytes[k] {
                    b'"' => val.push('"'),
                    b'\\' => val.push('\\'),
                    b'/' => val.push('/'),
                    b'n' => val.push('\n'),
                    b't' => val.push('\t'),
                    b'r' => val.push('\r'),
                    b'b' => val.push('\u{0008}'),
                    b'f' => val.push('\u{000C}'),
                    b'u' => {
                        // \uXXXX
                        let hex = &rest[k + 1..k + 5];
                        if let Ok(cp) = u32::from_str_radix(hex, 16) {
                            if let Some(c) = char::from_u32(cp) {
                                val.push(c);
                            }
                        }
                        k += 4;
                    }
                    other => val.push(other as char),
                }
            } else if ch == b'"' {
                return Some(val);
            } else {
                // Multi-byte UTF-8: copy the full char.
                let c = rest[k..].chars().next().unwrap();
                val.push(c);
                k += c.len_utf8();
                continue;
            }
            k += 1;
        }
        None
    }
}
