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
//!
//! Implementation note: this is a direct hand-roll of Biopython's `PDBParser._parse_coordinates`
//! and `MMCIFParser._build_structure` (DockQ's subclasses), NOT a `pdbtbx` transcode. pdbtbx
//! builds its own normalized model (its own altloc/element/model/chain-ordering rules); matching
//! Biopython byte-for-byte (auth ids, file-order chains, the exact `seq1` table, the `line[12:16]`
//! space rule, exact float32 from the literal coordinate text) is fewer moving parts and lower
//! risk done directly against the text. The Python source is mirrored line-for-line below; the
//! oracle JSON dumps (`oracle/dumps/parse{,_het}/`, regenerable via `oracle/dump_parse.py`) are
//! the judge — validated to match EXACTLY (chain order, sequence, residue ids, per-residue atom
//! names+order+count, and every coordinate to the f32 bit pattern) for all 13 example files on
//! both the core (`parse_hetatms=false`) and small-molecule (`parse_hetatms=true`) paths.

use std::io::Read;

use indexmap::IndexMap;

use crate::error::{DockQError, Result};
use crate::model::{Atom, Chain, Residue, Structure};

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
    let text = read_text(path)?;

    // Format detection by content (NO silent fallback). DockQ's `load_PDB` tries the
    // PDBParser first and falls back to MMCIFParser on *any* exception; we instead detect
    // explicitly. mmCIF is identified by the presence of `_atom_site.` loop tags; PDB by
    // ATOM/HETATM coordinate records. Neither => UnknownFormat.
    let is_mmcif = text.contains("_atom_site.");
    let has_pdb_atoms = text.lines().any(|l| {
        let r = record_type(l);
        r == "ATOM  " || r == "HETATM"
    });

    let mut structure = if is_mmcif {
        parse_mmcif(path, &text, chains, parse_hetatms, model_number)?
    } else if has_pdb_atoms {
        parse_pdb(path, &text, chains, parse_hetatms, model_number)?
    } else {
        return Err(DockQError::UnknownFormat {
            path: path.to_string(),
            detail: "no `_atom_site.` mmCIF tags and no ATOM/HETATM records found".to_string(),
        });
    };
    structure.id = path.to_string();
    Ok(structure)
}

// ---------------------------------------------------------------------------
// I/O + gzip
// ---------------------------------------------------------------------------

fn read_text(path: &str) -> Result<String> {
    let bytes = std::fs::read(path).map_err(|source| DockQError::Io {
        path: path.to_string(),
        source,
    })?;
    let gz = path.ends_with(".gz") || (bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b);
    let raw = if gz {
        let mut d = flate2::read::GzDecoder::new(&bytes[..]);
        let mut s = Vec::new();
        d.read_to_end(&mut s).map_err(|source| DockQError::Io {
            path: path.to_string(),
            source,
        })?;
        s
    } else {
        bytes
    };
    // Biopython reads in text mode; structure files are ASCII/UTF-8.
    Ok(String::from_utf8_lossy(&raw).into_owned())
}

// ---------------------------------------------------------------------------
// seq1 (Bio.SeqUtils.seq1) with DockQ's custom_map = {MSE:M, CME:C}.
// Resolved table = protein_letters_3to1_extended (upper-cased) overridden by custom_map.
// Unknown 3-letter -> "X". (Extracted from Biopython 1.87 in the baseline venv.)
// ---------------------------------------------------------------------------

fn seq1_three(resname_upper: &str) -> &'static str {
    match resname_upper {
        "ALA" => "A",
        "ARG" => "R",
        "ASN" => "N",
        "ASP" => "D",
        "ASX" => "B",
        "CME" => "C",
        "CYS" => "C",
        "GLN" => "Q",
        "GLU" => "E",
        "GLX" => "Z",
        "GLY" => "G",
        "HIS" => "H",
        "ILE" => "I",
        "LEU" => "L",
        "LYS" => "K",
        "MET" => "M",
        "MSE" => "M",
        "PHE" => "F",
        "PRO" => "P",
        "PYL" => "O",
        "SEC" => "U",
        "SER" => "S",
        "THR" => "T",
        "TRP" => "W",
        "TYR" => "Y",
        "VAL" => "V",
        "XAA" => "X",
        "XLE" => "J",
        _ => "X",
    }
}

/// Capitalized element symbols present in `Bio.Data.IUPACData.atom_weights` (110 keys).
/// Used by `assign_element` to validate a guessed element exactly as Biopython does.
const ATOM_WEIGHT_ELEMENTS: &[&str] = &[
    "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "B", "Ba", "Be", "Bh", "Bi", "Bk", "Br", "C",
    "Ca", "Cd", "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "D", "Db", "Dy", "Er", "Es", "Eu",
    "F", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge", "H", "He", "Hf", "Hg", "Ho", "Hs", "I", "In", "Ir",
    "K", "Kr", "La", "Li", "Lr", "Lu", "Md", "Mg", "Mn", "Mo", "Mt", "N", "Na", "Nb", "Nd", "Ne",
    "Ni", "No", "Np", "O", "Os", "P", "Pa", "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra", "Rb",
    "Re", "Rf", "Rh", "Rn", "Ru", "S", "Sb", "Sc", "Se", "Sg", "Si", "Sm", "Sn", "Sr", "Ta", "Tb",
    "Tc", "Te", "Th", "Ti", "Tl", "Tm", "U", "V", "W", "Xe", "Y", "Yb", "Zn", "Zr",
];

/// Python `str.capitalize()` for ASCII: first char upper, rest lower.
fn capitalize_ascii(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for (i, c) in s.chars().enumerate() {
        if i == 0 {
            out.extend(c.to_uppercase());
        } else {
            out.extend(c.to_lowercase());
        }
    }
    out
}

fn is_known_element(putative: &str) -> bool {
    let cap = capitalize_ascii(putative);
    ATOM_WEIGHT_ELEMENTS.iter().any(|&e| e == cap)
}

/// Replicates `Bio.PDB.Atom.Atom._assign_element`. `passed` is the (already stripped+upper-cased)
/// element column; `name` is the space-handled atom name; `fullname` is the raw `line[12:16]`.
/// Returns the element that Biopython would store on the Atom (guessed from the name when the
/// passed value is blank or not a recognized element). Not re-cased: returns the putative token
/// as-is, or "X" when unidentifiable — matching Biopython exactly.
fn assign_element(passed: &str, name: &str, fullname: &str) -> String {
    let needs_guess = passed.is_empty() || !is_known_element(passed);
    if !needs_guess {
        return passed.to_string();
    }
    let fb = fullname.as_bytes();
    let full0_alpha = fb.first().map(|c| c.is_ascii_alphabetic()).unwrap_or(false);
    // fullname[2:].isdigit(): Python's str.isdigit() is False for the empty string.
    let full2 = if fullname.len() > 2 { &fullname[2..] } else { "" };
    let full2_isdigit = !full2.is_empty() && full2.chars().all(|c| c.is_ascii_digit());

    let name_trim = name.trim();
    let nb = name.as_bytes();
    let putative: String = if full0_alpha && !full2_isdigit {
        name_trim.to_string()
    } else if nb.first().map(|c| c.is_ascii_digit()).unwrap_or(false) {
        // name[1]
        name.chars().nth(1).map(|c| c.to_string()).unwrap_or_default()
    } else {
        // name[0]
        name.chars().next().map(|c| c.to_string()).unwrap_or_default()
    };

    if is_known_element(&putative) {
        putative
    } else {
        "X".to_string()
    }
}

/// DockQ's per-residue one-letter contribution (only for standard / non-het residues):
///   len==3 -> seq1(resname, custom_map); len==2 -> resname[-1]; else -> resname.
fn residue_to_one(resname: &str) -> String {
    let n = resname.chars().count();
    if n == 3 {
        seq1_three(&resname.to_ascii_uppercase()).to_string()
    } else if n == 2 {
        // resname[-1]: last character, verbatim (no upper-casing).
        resname.chars().last().unwrap().to_string()
    } else {
        resname.to_string()
    }
}

// ---------------------------------------------------------------------------
// StructureBuilder: replicates the subset of Biopython's StructureBuilder semantics
// that DockQ observes through `get_atoms()` / `get_unpacked_list()` / chain.sequence.
//
// Per-residue atoms are deduplicated to ONE per atom-name:
//   * insertion order preserved (the order of *first* appearance of each name);
//   * the representative for a name is the highest-occupancy altloc, ties broken by
//     first-in-file (Biopython `DisorderedAtom.disordered_add`: replace only when
//     `occupancy > last_occupancy`, strictly greater).
// This makes `atoms.len()` == Biopython `len(list(res.get_atoms()))`
//   == `len(set(a.id for a in res.get_unpacked_list()))` (both dedup by name).
// ---------------------------------------------------------------------------

struct BResidue {
    het_flag: char,
    resseq: i64,
    icode: char,
    resname: String,
    resname1: String,
    /// (atom, occupancy_of_selected_child) per slot, in insertion order.
    atoms: Vec<(Atom, f64)>,
    /// name -> index into `atoms`.
    atom_index: IndexMap<String, usize>,
    /// Residue identity key used to detect residue changes: (het_flag, resseq, icode).
    key: (char, i64, char),
}

struct BChain {
    id: String,
    residues: Vec<BResidue>,
    /// One-letter sequence, built incrementally as residues are opened (mirrors DockQ's
    /// `sequences[chain]`): standard residues append their `seq1` code; a het residue sets it
    /// to the het resname.
    sequence: String,
    is_het: Option<String>,
}

#[derive(Default)]
struct BModel {
    chains: Vec<BChain>,
    chain_index: IndexMap<String, usize>,
}

#[derive(Default)]
struct Builder {
    models: Vec<BModel>,
    cur_model: Option<usize>,
    cur_chain: Option<usize>,
    cur_residue: Option<usize>,
}

impl Builder {
    fn init_model(&mut self) {
        self.models.push(BModel::default());
        self.cur_model = Some(self.models.len() - 1);
        self.cur_chain = None;
        self.cur_residue = None;
    }

    fn model(&mut self) -> &mut BModel {
        let i = self.cur_model.expect("model open");
        &mut self.models[i]
    }

    /// init_chain: reuse existing chain object with same id (discontinuous chain), else create.
    fn init_chain(&mut self, id: &str) {
        let m = self.model();
        if let Some(&idx) = m.chain_index.get(id) {
            self.cur_chain = Some(idx);
        } else {
            m.chains.push(BChain {
                id: id.to_string(),
                residues: Vec::new(),
                sequence: String::new(),
                is_het: None,
            });
            let idx = m.chains.len() - 1;
            m.chain_index.insert(id.to_string(), idx);
            self.cur_chain = Some(idx);
        }
        self.cur_residue = None;
    }

    fn chain(&mut self) -> &mut BChain {
        let mi = self.cur_model.expect("model open");
        let ci = self.cur_chain.expect("chain open");
        &mut self.models[mi].chains[ci]
    }

    /// init_residue: if the current residue already matches this key+resname, keep it;
    /// otherwise open a new residue object (append). Mirrors Biopython appending a new
    /// Residue on every (id,resname) change. (Full DisorderedResidue point-mutation
    /// handling is unused by the example oracle set; a recurring key creates a new
    /// residue object here, which matches the dump's residue iteration count.)
    fn init_residue(&mut self, resname: &str, het_flag: char, resseq: i64, icode: char) {
        let key = (het_flag, resseq, icode);
        let ch = self.chain();
        let new_res = BResidue {
            het_flag,
            resseq,
            icode,
            resname: resname.to_string(),
            resname1: String::new(),
            atoms: Vec::new(),
            atom_index: IndexMap::new(),
            key,
        };
        ch.residues.push(new_res);
        let idx = ch.residues.len() - 1;
        self.cur_residue = Some(idx);
    }

    fn residue(&mut self) -> &mut BResidue {
        let mi = self.cur_model.expect("model open");
        let ci = self.cur_chain.expect("chain open");
        let ri = self.cur_residue.expect("residue open");
        &mut self.models[mi].chains[ci].residues[ri]
    }

    /// init_atom with disordered-atom dedup. `altloc == ' '` => ordered atom; a duplicate
    /// ordered name is a Biopython "defined twice" error which DockQ (QUIET) swallows and
    /// drops — we keep the first and ignore the later one. `altloc != ' '` => disordered;
    /// representative is highest occupancy, ties first-in-file.
    fn init_atom(&mut self, name: &str, element: &str, altloc: char, coord: [f32; 3], occupancy: f64) {
        let atom = Atom {
            name: name.to_string(),
            element: element.to_string(),
            altloc,
            coord,
        };
        let res = self.residue();
        if let Some(&slot) = res.atom_index.get(name) {
            // Name already present in this residue.
            let (existing, last_occ) = &mut res.atoms[slot];
            if altloc != ' ' {
                // Disordered: replace selected child iff occupancy strictly greater.
                if occupancy > *last_occ {
                    *existing = atom;
                    *last_occ = occupancy;
                }
            }
            // else: ordered duplicate name -> Biopython raises & QUIET drops; keep first.
        } else {
            let slot = res.atoms.len();
            res.atoms.push((atom, occupancy));
            res.atom_index.insert(name.to_string(), slot);
        }
    }

    fn finish(self, model_number: usize) -> Result<IndexMap<String, Chain>> {
        let available = self.models.len();
        let m = self
            .models
            .into_iter()
            .nth(model_number)
            .ok_or(DockQError::ModelOutOfRange {
                requested: model_number,
                available,
            })?;
        let mut out: IndexMap<String, Chain> = IndexMap::new();
        for bc in m.chains {
            let residues = bc
                .residues
                .into_iter()
                .map(|br| Residue {
                    het_flag: br.het_flag,
                    resseq: br.resseq,
                    icode: br.icode,
                    resname: br.resname,
                    resname1: br.resname1,
                    atoms: br.atoms.into_iter().map(|(a, _)| a).collect(),
                })
                .collect();
            out.insert(
                bc.id.clone(),
                Chain {
                    id: bc.id,
                    residues,
                    sequence: bc.sequence,
                    is_het: bc.is_het,
                },
            );
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------
// PDB parser — mirrors DockQ.PDBParser._parse_coordinates
// ---------------------------------------------------------------------------

#[inline]
fn record_type(line: &str) -> String {
    // Python `line[0:6]`: first 6 chars, right-padded conceptually. We compare against
    // 6-char literals, so take exactly the first 6 bytes (ASCII) without padding;
    // shorter lines simply won't match "ATOM  "/"HETATM".
    let b = line.as_bytes();
    let end = b.len().min(6);
    // Pad to 6 with spaces to match Python slice semantics for short lines used elsewhere.
    let mut s = String::with_capacity(6);
    s.push_str(std::str::from_utf8(&b[..end]).unwrap_or(""));
    while s.len() < 6 {
        s.push(' ');
    }
    s
}

/// Python `line[a:b]` on a (already newline-stripped) line: byte slice, clamped, never panics.
#[inline]
fn slice(line: &str, a: usize, b: usize) -> &str {
    let bytes = line.as_bytes();
    let lo = a.min(bytes.len());
    let hi = b.min(bytes.len());
    if lo >= hi {
        ""
    } else {
        std::str::from_utf8(&bytes[lo..hi]).unwrap_or("")
    }
}

/// Python `line[i]` single char (as &str); empty if out of range.
#[inline]
fn at(line: &str, i: usize) -> char {
    line.as_bytes().get(i).map(|&c| c as char).unwrap_or(' ')
}

fn parse_pdb(
    path: &str,
    text: &str,
    chains: &[String],
    parse_hetatms: bool,
    model_number: usize,
) -> Result<Structure> {
    let mut b = Builder::default();
    let mut model_open = false;

    for raw_line in text.split_inclusive('\n') {
        // Python iterates readlines() then does line.rstrip("\n"). We also drop a trailing
        // '\r' so CRLF files behave like Python text-mode (which strips \r on read).
        let line = raw_line.trim_end_matches('\n').trim_end_matches('\r');
        let rt = record_type(line);

        if line.trim().is_empty() {
            continue; // skip empty lines
        } else if rt == "HETATM" && !parse_hetatms {
            continue;
        } else if rt == "ATOM  " || rt == "HETATM" {
            if !model_open {
                b.init_model();
                model_open = true;
            }
            let chainid = at(line, 21);
            let chainid_s = chainid.to_string();
            if !chains.is_empty() && !chains.iter().any(|c| c == &chainid_s) {
                continue;
            }
            // Column element, stripped+upper. Used ONLY for the hydrogen skip-check (Biopython
            // skips before constructing the Atom, so the element *guess* never reinstates an H).
            let col_element = slice(line, 76, 78).trim().to_ascii_uppercase();
            if col_element == "H" {
                continue;
            }
            // fullname = line[12:16]; if it has internal spaces (split != 1 token) keep as-is,
            // else the single stripped token.
            let fullname = slice(line, 12, 16);
            let split: Vec<&str> = fullname.split_whitespace().collect();
            let name: String = if split.len() != 1 {
                fullname.to_string()
            } else {
                split[0].to_string()
            };
            // Stored element: Biopython's Atom._assign_element guesses from the name when the
            // column element is blank or unrecognized (e.g. coordinate-only files like model.pdb).
            let element = assign_element(&col_element, &name, fullname);
            let altloc = at(line, 16);
            let resname = slice(line, 17, 20).trim().to_string();
            let het_flag = if rt == "HETATM" { 'H' } else { ' ' };

            // resseq = int(line[22:26].split()[0]); icode = line[26].
            let resseq_field = slice(line, 22, 26);
            let resseq_tok = match resseq_field.split_whitespace().next() {
                Some(t) => t,
                None => {
                    return Err(DockQError::Parse {
                        path: path.to_string(),
                        msg: format!("missing residue sequence number in line: {line:?}"),
                    })
                }
            };
            let resseq: i64 = resseq_tok.parse().map_err(|_| DockQError::Parse {
                path: path.to_string(),
                msg: format!("invalid residue sequence number {resseq_tok:?}"),
            })?;
            let icode = at(line, 26);

            // coordinates: float(text) then cast to f32 (Biopython np.array(...,'f')).
            let x = parse_coord(path, slice(line, 30, 38))?;
            let y = parse_coord(path, slice(line, 38, 46))?;
            let z = parse_coord(path, slice(line, 46, 54))?;

            // occupancy: float(line[54:60]); missing -> None in Biopython (we use 0.0 only
            // as the altloc tiebreak weight; the example set has no altlocs).
            let occupancy = slice(line, 54, 60).trim().parse::<f64>().unwrap_or(0.0);

            // chain / residue bookkeeping (mirrors the init_chain / init_residue logic).
            let new_chain = b.cur_chain.is_none() || b.chain().id != chainid_s;
            if new_chain {
                b.init_chain(&chainid_s);
                b.init_residue(&resname, het_flag, resseq, icode);
                set_residue_sequence(&mut b, het_flag, &resname, /*append=*/ false);
            } else {
                let key = (het_flag, resseq, icode);
                let changed = {
                    let r = b.residue();
                    r.key != key || r.resname != resname
                };
                if changed {
                    b.init_residue(&resname, het_flag, resseq, icode);
                    set_residue_sequence(&mut b, het_flag, &resname, /*append=*/ true);
                }
            }

            b.init_atom(&name, &element, altloc, [x, y, z], occupancy);
        } else if rt == "MODEL " {
            // New explicit model.
            b.init_model();
            model_open = true;
        } else if rt == "END   " || rt == "CONECT" {
            break; // end of atomic data
        } else if rt == "ENDMDL" {
            model_open = false;
            b.cur_chain = None;
            b.cur_residue = None;
        }
        // ANISOU / SIGUIJ / SIGATM / TER / unrecognized: ignored (no atoms produced).
    }

    let chains_map = b.finish(model_number)?;
    Ok(Structure {
        chains: chains_map,
        id: path.to_string(),
    })
}

fn parse_coord(path: &str, field: &str) -> Result<f32> {
    let t = field.trim();
    let v: f64 = t.parse().map_err(|_| DockQError::Parse {
        path: path.to_string(),
        msg: format!("invalid coordinate {t:?}"),
    })?;
    Ok(v as f32)
}

/// Replicates the sequence/is_het update DockQ performs at residue init:
///   het_flag == ' ': append (or set, if first in chain) the one-letter code;
///   het_flag != ' ': sequence = resname; is_het = resname.
fn set_residue_sequence(b: &mut Builder, het_flag: char, resname: &str, append: bool) {
    if het_flag == ' ' {
        let one = residue_to_one(resname);
        {
            let r = b.residue();
            r.resname1 = one.clone();
        }
        let ch = b.chain();
        if append {
            ch.sequence.push_str(&one);
        } else {
            ch.sequence = one;
        }
    } else {
        let resname_owned = resname.to_string();
        {
            let r = b.residue();
            r.resname1 = String::new();
        }
        let ch = b.chain();
        ch.sequence = resname_owned.clone();
        ch.is_het = Some(resname_owned);
    }
}

// ---------------------------------------------------------------------------
// mmCIF parser — mirrors DockQ.MMCIFParser._build_structure (auth chains + auth residues)
// ---------------------------------------------------------------------------

/// Column indices into an `_atom_site` row for the fields DockQ reads.
struct AtomSiteCols {
    group_pdb: usize,
    type_symbol: Option<usize>,
    label_atom_id: usize,
    label_comp_id: usize,
    auth_asym_id: usize,  // chain (auth_chains=True)
    label_asym_id: usize, // fallback chain
    cartn_x: usize,
    cartn_y: usize,
    cartn_z: usize,
    label_alt_id: usize,
    auth_seq_id: Option<usize>, // resseq (auth_residues=True), falls back to label_seq_id
    label_seq_id: usize,
    pdbx_ins_code: usize,
    occupancy: usize,
    model_num: Option<usize>,
}

fn parse_mmcif(
    path: &str,
    text: &str,
    chains: &[String],
    parse_hetatms: bool,
    model_number: usize,
) -> Result<Structure> {
    let (cols, rows) = extract_atom_site_loop(path, text)?;

    // DockQ's `load_PDB` passes `auth_chains = not small_molecule`, i.e. `auth_chains == !parse_hetatms`:
    // the core path (parse_hetatms=false) keys chains by auth_asym_id; the small-molecule path
    // (parse_hetatms=true) keys chains by label_asym_id (so het groups like NDP/HOH become their
    // own label chains C/D/E rather than collapsing into the polymer's auth chain). `auth_residues`
    // is always True in DockQ's subclass, so resseq uses auth_seq_id (falling back to label_seq_id).
    let use_auth_chain = !parse_hetatms;
    let use_auth_residue = cols.auth_seq_id.is_some(); // auth_residues=True; fall back to label

    let mut b = Builder::default();
    let mut current_serial: Option<String> = None;
    let mut model_initialized = false;

    for row in &rows {
        let get = |i: usize| -> &str { row.get(i).map(|s| s.as_str()).unwrap_or("") };

        let chainid = if use_auth_chain {
            get(cols.auth_asym_id)
        } else {
            get(cols.label_asym_id)
        };
        if !chains.is_empty() && !chains.iter().any(|c| c == chainid) {
            continue;
        }
        let fieldname = get(cols.group_pdb);
        if fieldname == "HETATM" && !parse_hetatms {
            continue;
        }
        // type_symbol (upper) drives the hydrogen skip-check; if absent Biopython passes None
        // (never "H", so nothing is skipped). Stored element is guessed from the name when
        // blank/unrecognized, exactly as the PDB path (Biopython mmCIF passes name as fullname).
        let col_element = cols
            .type_symbol
            .map(|i| get(i).to_ascii_uppercase())
            .unwrap_or_default();
        if cols.type_symbol.is_some() && col_element == "H" {
            continue;
        }

        let resname = get(cols.label_comp_id).to_string();
        let mut altloc = get(cols.label_alt_id).chars().next().unwrap_or(' ');
        if altloc == '.' || altloc == '?' {
            altloc = ' ';
        }
        let resseq_raw = if use_auth_residue {
            get(cols.auth_seq_id.unwrap())
        } else {
            get(cols.label_seq_id)
        };
        if resseq_raw == "." {
            // Non-existing residue ID -> Biopython warns & continues (skip atom).
            continue;
        }
        let int_resseq: i64 = resseq_raw.parse().map_err(|_| DockQError::Parse {
            path: path.to_string(),
            msg: format!("invalid mmCIF residue id {resseq_raw:?}"),
        })?;
        let mut icode = get(cols.pdbx_ins_code).chars().next().unwrap_or(' ');
        if icode == '.' || icode == '?' {
            icode = ' ';
        }
        let name = get(cols.label_atom_id).to_string();
        // Biopython mmCIF passes `name` as the fullname to Atom(); _assign_element guesses from
        // it when type_symbol is blank/unrecognized (col_element is "" when the column is absent).
        let element = assign_element(&col_element, &name, &name);

        let x = parse_coord(path, get(cols.cartn_x))?;
        let y = parse_coord(path, get(cols.cartn_y))?;
        let z = parse_coord(path, get(cols.cartn_z))?;
        let occupancy = get(cols.occupancy).trim().parse::<f64>().unwrap_or(0.0);

        let het_flag = if fieldname == "HETATM" { 'H' } else { ' ' };

        // Model handling: if a model column exists, a change in its value starts a new model
        // (Biopython increments the array index). Otherwise a single model.
        if let Some(mi) = cols.model_num {
            let serial = get(mi).to_string();
            if current_serial.as_deref() != Some(serial.as_str()) {
                current_serial = Some(serial);
                b.init_model();
                model_initialized = true;
            }
        } else if !model_initialized {
            b.init_model();
            model_initialized = true;
        }

        // chain bookkeeping
        let new_chain = b.cur_chain.is_none() || b.chain().id != chainid;
        if new_chain {
            b.init_chain(chainid);
        }

        // residue bookkeeping: (current_residue_id != resseq) OR (current_resname != resname)
        let key = (het_flag, int_resseq, icode);
        let need_new_res = match b.cur_residue {
            None => true,
            Some(_) => {
                let r = b.residue();
                r.key != key || r.resname != resname
            }
        };
        if need_new_res {
            b.init_residue(&resname, het_flag, int_resseq, icode);
            // mmCIF: Biopython initializes the chain's sequence to "" when the chain is first
            // seen, then always `+=`'s for standard residues -> always append here.
            set_residue_sequence(&mut b, het_flag, &resname, /*append=*/ true);
        }

        b.init_atom(&name, &element, altloc, [x, y, z], occupancy);
    }

    if !model_initialized {
        return Err(DockQError::Parse {
            path: path.to_string(),
            msg: "mmCIF `_atom_site` loop produced no atoms".to_string(),
        });
    }

    let chains_map = b.finish(model_number)?;
    Ok(Structure {
        chains: chains_map,
        id: path.to_string(),
    })
}

/// Extract the `_atom_site` loop: returns the resolved column map and the data rows
/// (each row tokenized, quotes handled). Mirrors what MMCIF2Dict yields for `_atom_site.*`.
fn extract_atom_site_loop(
    path: &str,
    text: &str,
) -> Result<(AtomSiteCols, Vec<Vec<String>>)> {
    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;
    // Find a `loop_` whose first tag line is `_atom_site.`.
    let mut header_tags: Vec<String> = Vec::new();
    let mut data_start: Option<usize> = None;
    while i < lines.len() {
        if lines[i].trim() == "loop_" {
            // collect following tag lines
            let mut j = i + 1;
            let mut tags: Vec<String> = Vec::new();
            while j < lines.len() {
                let t = lines[j].trim();
                if t.starts_with('_') {
                    tags.push(t.to_string());
                    j += 1;
                } else {
                    break;
                }
            }
            if tags.first().map(|t| t.starts_with("_atom_site.")).unwrap_or(false) {
                header_tags = tags;
                data_start = Some(j);
                break;
            }
            i = j;
        } else {
            i += 1;
        }
    }

    let data_start = match data_start {
        Some(d) => d,
        None => {
            return Err(DockQError::Parse {
                path: path.to_string(),
                msg: "no `_atom_site` loop found in mmCIF".to_string(),
            })
        }
    };

    // Map tag -> column index.
    let idx_of = |needle: &str| -> Option<usize> {
        header_tags.iter().position(|t| t == needle)
    };
    let req = |needle: &str| -> Result<usize> {
        idx_of(needle).ok_or_else(|| DockQError::Parse {
            path: path.to_string(),
            msg: format!("mmCIF `_atom_site` loop missing column {needle}"),
        })
    };

    let cols = AtomSiteCols {
        group_pdb: req("_atom_site.group_PDB")?,
        type_symbol: idx_of("_atom_site.type_symbol"),
        label_atom_id: req("_atom_site.label_atom_id")?,
        label_comp_id: req("_atom_site.label_comp_id")?,
        auth_asym_id: req("_atom_site.auth_asym_id")?,
        label_asym_id: req("_atom_site.label_asym_id")?,
        cartn_x: req("_atom_site.Cartn_x")?,
        cartn_y: req("_atom_site.Cartn_y")?,
        cartn_z: req("_atom_site.Cartn_z")?,
        label_alt_id: req("_atom_site.label_alt_id")?,
        auth_seq_id: idx_of("_atom_site.auth_seq_id"),
        label_seq_id: req("_atom_site.label_seq_id")?,
        pdbx_ins_code: req("_atom_site.pdbx_PDB_ins_code")?,
        occupancy: req("_atom_site.occupancy")?,
        model_num: idx_of("_atom_site.pdbx_PDB_model_num"),
    };

    let ncols = header_tags.len();
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut k = data_start;
    while k < lines.len() {
        let line = lines[k];
        let trimmed = line.trim_start();
        // End of the atom_site data block: next loop_, next category, comment, or `data_`/`save_`.
        if trimmed == "loop_"
            || trimmed.starts_with('_')
            || trimmed.starts_with('#')
            || trimmed.starts_with("data_")
            || trimmed.starts_with("save_")
        {
            break;
        }
        if trimmed.is_empty() {
            k += 1;
            continue;
        }
        // Multiline `;` text fields do not occur inside `_atom_site` rows for these files;
        // if encountered we'd mis-tokenize, so guard against it explicitly.
        if trimmed.starts_with(';') {
            return Err(DockQError::Parse {
                path: path.to_string(),
                msg: "unexpected multiline `;` field inside `_atom_site` loop".to_string(),
            });
        }
        let toks = tokenize_cif_line(line);
        if toks.len() == ncols {
            rows.push(toks);
        } else if !toks.is_empty() {
            // A data row whose token count doesn't match the column count means either a
            // wrapped row (values spanning multiple lines) or corruption. The example files
            // are one-row-per-line; surface the mismatch rather than silently misalign.
            return Err(DockQError::Parse {
                path: path.to_string(),
                msg: format!(
                    "mmCIF `_atom_site` row has {} tokens, expected {}: {:?}",
                    toks.len(),
                    ncols,
                    line
                ),
            });
        }
        k += 1;
    }

    Ok((cols, rows))
}

/// Tokenize one mmCIF data line into values, honoring single/double quoting (whitespace
/// inside quotes is preserved; the quotes themselves are stripped, matching MMCIF2Dict).
fn tokenize_cif_line(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let bytes = line.as_bytes();
    let mut i = 0;
    let n = bytes.len();
    while i < n {
        // skip whitespace
        while i < n && (bytes[i] == b' ' || bytes[i] == b'\t' || bytes[i] == b'\r') {
            i += 1;
        }
        if i >= n {
            break;
        }
        let c = bytes[i];
        if c == b'\'' || c == b'"' {
            // Quoted value: ends at the matching quote *followed by whitespace or EOL*
            // (mmCIF rule). We use the simpler common case: quote ... quote.
            let quote = c;
            i += 1;
            let start = i;
            while i < n {
                if bytes[i] == quote && (i + 1 >= n || bytes[i + 1] == b' ' || bytes[i + 1] == b'\t')
                {
                    break;
                }
                i += 1;
            }
            out.push(String::from_utf8_lossy(&bytes[start..i.min(n)]).into_owned());
            if i < n {
                i += 1; // skip closing quote
            }
        } else {
            let start = i;
            while i < n && bytes[i] != b' ' && bytes[i] != b'\t' && bytes[i] != b'\r' {
                i += 1;
            }
            out.push(String::from_utf8_lossy(&bytes[start..i]).into_owned());
        }
    }
    out
}

// ===========================================================================
// Tests — assert key invariants against hardcoded oracle-derived values.
// (Expected values cross-checked against oracle/dumps/parse/*.json.)
// ===========================================================================
#[cfg(test)]
mod tests {
    use super::*;

    // Vendored example structures (a small MIT-licensed subset of upstream DockQ's
    // examples/), resolved at compile time so `cargo test` passes on any machine.
    const EX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data");

    fn load(name: &str) -> Structure {
        load_structure(&format!("{EX}/{name}"), &[], false, 0).expect("load ok")
    }

    #[test]
    fn pdb_1a2k_chains_and_residues() {
        let s = load("1A2K_r_l_b.pdb");
        // chain ids + order
        assert_eq!(s.chain_ids(), vec!["A", "B", "C"]);
        // per-chain residue counts (from oracle dump)
        assert_eq!(s.chain("A").unwrap().residues.len(), 124);
        assert_eq!(s.chain("B").unwrap().residues.len(), 124);
        assert_eq!(s.chain("C").unwrap().residues.len(), 196);
        // sequence length == residue count (polymer invariant)
        for id in ["A", "B", "C"] {
            let c = s.chain(id).unwrap();
            assert_eq!(c.sequence.chars().count(), c.residues.len(), "chain {id}");
        }
        // first residue of A = LYS 4, 9 atoms; first atom N @ exact f32 bits
        let a0 = &s.chain("A").unwrap().residues[0];
        assert_eq!(a0.resname, "LYS");
        assert_eq!(a0.resseq, 4);
        assert_eq!(a0.icode, ' ');
        assert_eq!(a0.het_flag, ' ');
        assert_eq!(a0.atoms.len(), 9);
        let n = &a0.atoms[0];
        assert_eq!(n.name, "N");
        assert_eq!(n.element, "N");
        // 28.189 / 5.020 / 62.680 as f32 (exact bit patterns).
        assert_eq!(n.coord[0].to_bits(), 28.189_f32.to_bits());
        assert_eq!(n.coord[1].to_bits(), 5.02_f32.to_bits());
        assert_eq!(n.coord[2].to_bits(), 62.68_f32.to_bits());
        // last residue of A = GLY 127
        let alast = s.chain("A").unwrap().residues.last().unwrap();
        assert_eq!(alast.resname, "GLY");
        assert_eq!(alast.resseq, 127);
    }

    #[test]
    fn pdb_1a2k_sequence_first_residues() {
        let s = load("1A2K_r_l_b.pdb");
        // chain A starts LYS(K) -> sequence begins with 'K'
        assert!(s.chain("A").unwrap().sequence.starts_with('K'));
        // chain C starts GLN(Q)
        assert!(s.chain("C").unwrap().sequence.starts_with('Q'));
    }

    #[test]
    fn hetatm_excluded_in_core_path() {
        // 1HHO has HEM (HETATM); core path (parse_hetatms=false) must exclude all het.
        let s = load("1HHO_hem.cif");
        for (_id, c) in &s.chains {
            assert!(c.is_het.is_none(), "no het chains in core path");
            for r in &c.residues {
                assert_eq!(r.het_flag, ' ', "no het residues in core path");
            }
            assert_eq!(c.sequence.chars().count(), c.residues.len());
        }
    }

    #[test]
    fn mmcif_1exb_gz_exact() {
        // 1EXB.cif.gz: gzip + mmCIF; auth_asym_id chains in file order = [A, E].
        let s = load_structure(&format!("{EX}/1EXB.cif.gz"), &[], false, 0).expect("load");
        assert_eq!(s.chain_ids(), vec!["A", "E"]);
        assert_eq!(s.chain("A").unwrap().residues.len(), 326);
        assert_eq!(s.chain("E").unwrap().residues.len(), 91);
        for id in ["A", "E"] {
            let c = s.chain(id).unwrap();
            assert_eq!(c.sequence.chars().count(), c.residues.len(), "chain {id}");
        }
        // chain A sequence prefix and first residue (auth_seq_id 36, LEU).
        assert!(s.chain("A").unwrap().sequence.starts_with("LQFYRNLGKS"));
        let r0 = &s.chain("A").unwrap().residues[0];
        assert_eq!(r0.resname, "LEU");
        assert_eq!(r0.resseq, 36);
        assert_eq!(r0.atoms.len(), 8);
        let n = &r0.atoms[0];
        assert_eq!(n.name, "N");
        assert_eq!(n.element, "N"); // from type_symbol
        assert_eq!(n.coord[0].to_bits(), 0xbf8dd2f2);
    }

    #[test]
    fn pdb_model_pdb_element_guessed_from_name() {
        // model.pdb has truncated lines (no element column); Biopython guesses element
        // from the atom name. chains = {A, B}; chain A first residue GLY 2, atom0 N guessed "N".
        let s = load("model.pdb");
        assert_eq!(s.chain_ids(), vec!["A", "B"]);
        let a = s.chain("A").unwrap();
        assert_eq!(a.residues.len(), 374);
        assert_eq!(a.sequence.chars().count(), 374);
        let r0 = &a.residues[0];
        assert_eq!(r0.resname, "GLY");
        assert_eq!(r0.resseq, 2);
        assert_eq!(r0.atoms.len(), 4);
        let n = &r0.atoms[0];
        assert_eq!(n.name, "N");
        assert_eq!(n.element, "N", "element guessed from atom name");
        assert_eq!(n.coord[0].to_bits(), 0x41d9020c);
        assert_eq!(n.coord[1].to_bits(), 0x41b66873);
        assert_eq!(n.coord[2].to_bits(), 0x4244a1cb);
    }

    #[test]
    fn assign_element_rules() {
        // Recognized column element passes through unchanged.
        assert_eq!(assign_element("N", "N", " N  "), "N");
        assert_eq!(assign_element("SE", "SE", "SE  "), "SE");
        // Blank column, " CA " -> fullname[0]==' ' (not alpha) -> name[0]='C'.
        assert_eq!(assign_element("", "CA", " CA "), "C");
        // Blank column, " N  " -> name[0]='N'.
        assert_eq!(assign_element("", "N", " N  "), "N");
        // 4-char name HE21: fullname[0]='H' alpha, fullname[2:]="21" all-digit -> else branch;
        // name[0]='H' not digit -> name[0]='H'. (Not skipped: column element was blank.)
        assert_eq!(assign_element("", "HE21", "HE21"), "H");
        // digit-first hydrogen name "1HB": fullname[0]='1' not alpha -> else; name[0]='1' digit
        // -> name[1]='H'.
        assert_eq!(assign_element("", "1HB", "1HB "), "H");
    }

    #[test]
    fn chains_filter_restricts() {
        let only_a = load_structure(
            &format!("{EX}/1A2K_r_l_b.pdb"),
            &["A".to_string()],
            false,
            0,
        )
        .expect("load");
        assert_eq!(only_a.chain_ids(), vec!["A"]);
    }

    #[test]
    fn unknown_format_errors() {
        let dir = std::env::temp_dir();
        let p = dir.join("dockq_parser_bogus.txt");
        std::fs::write(&p, "this is not a structure file\n").unwrap();
        let r = load_structure(p.to_str().unwrap(), &[], false, 0);
        assert!(matches!(r, Err(DockQError::UnknownFormat { .. })));
        let _ = std::fs::remove_file(&p);
    }
}
