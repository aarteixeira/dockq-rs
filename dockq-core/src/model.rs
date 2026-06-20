//! Core data model. Mirrors the subset of Biopython's Structure/Chain/Residue/Atom
//! hierarchy that DockQ actually uses.
//!
//! Coordinates are **f32** on purpose: Biopython stores atom coords as numpy `'f'`
//! (float32) and the reference Cython `residue_distances` kernel computes in float32.
//! Matching that precision keeps us numerically aligned with the oracle (and is faster).

use indexmap::IndexMap;

/// A single heavy atom (hydrogens are dropped at parse time, matching DockQ).
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    /// Atom name as DockQ uses it (`atom.id`), e.g. "CA", "N", "OP1".
    /// Whitespace-stripped per the parser rules.
    pub name: String,
    /// Element symbol, uppercased (e.g. "C", "N", "SE"). Used to skip hydrogens.
    pub element: String,
    /// Alternate location indicator (' ' if none).
    pub altloc: char,
    /// Cartesian coordinates in Å, float32 (matches Biopython + the Cython kernel).
    pub coord: [f32; 3],
}

/// A residue (standard amino acid / nucleotide for the core path; het groups when
/// `parse_hetatms` is enabled).
#[derive(Clone, Debug, PartialEq)]
pub struct Residue {
    /// Hetero flag: ' ' for standard ATOM residues, 'H' for HETATM het groups.
    pub het_flag: char,
    /// Residue sequence number (auth_seq_id for mmCIF).
    pub resseq: i64,
    /// Insertion code (' ' if none).
    pub icode: char,
    /// Residue name, e.g. "ALA", "DA", "HEM".
    pub resname: String,
    /// One-letter code this residue contributed to the chain sequence (debug/reference).
    pub resname1: String,
    /// Heavy atoms, **deduplicated to one per atom name** (the Biopython `get_atoms`
    /// representative for altloc groups). File order is preserved. This guarantees
    /// `atoms.len()` equals the per-residue atom count fed to `residue_distances`.
    pub atoms: Vec<Atom>,
}

impl Residue {
    /// Number of (unique-named) atoms.
    #[inline]
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// First atom with the given name (the `subset_atoms` lookup for backbone atoms).
    pub fn atom_by_name(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.name == name)
    }
}

/// A polymer (or het) chain.
#[derive(Clone, Debug, PartialEq)]
pub struct Chain {
    pub id: String,
    pub residues: Vec<Residue>,
    /// One-letter sequence (concatenated `resname1` of standard residues), used for
    /// alignment. Invariant for polymer chains: `sequence.chars().count() == residues.len()`.
    pub sequence: String,
    /// Het identity: `None` for polymer chains; `Some(resname)` for het groups.
    pub is_het: Option<String>,
}

impl Chain {
    #[inline]
    pub fn len(&self) -> usize {
        self.residues.len()
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }
}

/// A parsed structure: one selected model, its chains in file order.
#[derive(Clone, Debug)]
pub struct Structure {
    /// Chains in file order, keyed by chain id.
    pub chains: IndexMap<String, Chain>,
    /// Source path / id (diagnostics).
    pub id: String,
}

impl Structure {
    pub fn chain(&self, id: &str) -> Option<&Chain> {
        self.chains.get(id)
    }
    pub fn chain_ids(&self) -> Vec<String> {
        self.chains.keys().cloned().collect()
    }
}

/// A formatted pairwise alignment, mirroring DockQ's `format_alignment` output.
#[derive(Clone, Debug, PartialEq)]
pub struct Alignment {
    /// Aligned model sequence (with '-' for gaps).
    pub seq_a: String,
    /// Match string: '|' identical column, '.' substitution, ' ' gap.
    pub matches: String,
    /// Aligned native sequence (with '-' for gaps).
    pub seq_b: String,
}

/// Dense row-major matrix of (squared) residue–residue distances (float32),
/// the output of `residue_distances`.
#[derive(Clone, Debug)]
pub struct DistMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f32>,
}

impl DistMatrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f32 {
        self.data[i * self.cols + j]
    }
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f32) {
        self.data[i * self.cols + j] = v;
    }
    /// Count of entries strictly below `threshold` (used for nat_total / clashes).
    pub fn count_below(&self, threshold: f32) -> u32 {
        self.data.iter().filter(|&&d| d < threshold).count() as u32
    }
}

/// Per-interface DockQ result (mirrors the `info` dict produced by `calc_DockQ`).
/// Scores are stored as f64 (geometry is computed in f32 then widened).
#[derive(Clone, Debug)]
pub struct InterfaceResult {
    pub dockq: f64,
    pub f1: f64,
    pub irmsd: f64,
    pub lrmsd: f64,
    pub fnat: f64,
    pub nat_correct: u32,
    pub nat_total: u32,
    pub fnonnat: f64,
    pub nonnat_count: u32,
    pub model_total: u32,
    pub clashes: u32,
    pub len1: usize,
    pub len2: usize,
    /// "receptor" / "ligand" for chain group 1.
    pub class1: String,
    /// "receptor" / "ligand" for chain group 2.
    pub class2: String,
    /// Model chain id mapped to native chain 1.
    pub chain1: String,
    /// Model chain id mapped to native chain 2.
    pub chain2: String,
}
