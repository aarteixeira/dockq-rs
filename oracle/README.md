# DockQ Ground-Truth Oracle

JSON dumps of the reference **DockQ v2.1.3** (Python + Cython) intermediate and
final values, for validating a Rust reimplementation (`dockq-rs`) byte-for-byte
or within tolerance.

The reference implementation is **correct by definition**. These dumps are
produced by calling the reference library's own functions (no DockQ logic is
reimplemented here), so they are authoritative.

- Reference repo: `/Users/Andre.Teixeira/projects/DockQ`
- Reference version: DockQ **2.1.3**, Python **3.12.8**, NumPy float32 coords
- Example structures: `/Users/Andre.Teixeira/projects/DockQ/examples/`
- Golden CLI outputs: `/Users/Andre.Teixeira/projects/DockQ/testdata/`

## Conventions (read first)

- **Distances are squared, in Angstrom^2.** Every residue-distance matrix DockQ
  produces holds the *squared* minimum inter-residue atom distance (the Cython
  kernel `residue_distances` never takes a sqrt). A 5 Angstrom contact is
  therefore `value < 25.0`; a 10 Angstrom interface is `value < 100.0`. All
  `n_lt_*` counts below are computed on the squared matrix against squared
  thresholds.
- **`coord` vs `coord_hex`.** Biopython stores atom coordinates as NumPy
  `float32`. `coord` is `[float(x), float(y), float(z)]` (the float32 value
  widened to a Python float for human reading — may show float64 rounding
  artifacts like `13.593000411987305`). `coord_hex` is the exact bit-faithful
  form: `[float.hex(float(c)) for c in atom.coord]`. To recover the exact
  float32 in Rust/Python: `np.float32(float.fromhex(h))`. Round-trip verified
  exact.
- **`bytes_sha256`.** SHA-256 of a residue-distance matrix serialized as
  *row-major (C-order) float32* via `matrix.astype(np.float32).tobytes()`. A
  Rust port that builds the same float32 matrix in row-major order and hashes
  its little-endian bytes will match this digest. Verified reproducible.
- **`het_flag`.** Biopython's residue hetero flag (`res.id[0]`). `" "` for
  standard ATOM residues; for HETATM it is the full form `"H_<RESNAME>"`
  (e.g. `"H_HEM"`, `"H_HOH"`).
- **Floats are full precision.** Nothing is rounded. Golden CLI prints 3
  decimals; these dumps keep all digits (e.g. iRMSD `1.4879852993911616e-06`).
- **All JSON is written with `sort_keys=True, indent=1`** (deterministic).

---

## (A) Parse dumps — `dumps/parse/` and `dumps/parse_het/`

One file per example structure, produced by
`DockQ.DockQ.load_PDB(path, small_molecule=False)` (`parse/`) and
`load_PDB(path, small_molecule=True)` (`parse_het/`).
Filename = `<original_filename>.json` (e.g. `1A2K_r_l_b.model.pdb.json`).
`small_molecule=True` additionally parses HETATM records, so `parse_het/`
dumps contain extra het chains (ligands, ions, waters) absent from `parse/`.

Errors are caught per-file and recorded in `"error"` (the run never aborts).

```jsonc
{
  "path": str,                  // absolute path passed to load_PDB
  "error": null | str,          // full Python traceback string, or null
  "chains": [
    {
      "id": str,                // chain.id (auth chain id)
      "sequence": str,          // chain.sequence; 1-letter for protein/NA,
                                //   or the het resname for a het chain
      "is_het": null | str,     // chain.is_het: null for polymer, else resname
                                //   (e.g. "HEM", "PO4", "HOH")
      "n_residues": int,        // number of residues iterated from the chain
      "residues": [
        {
          "resname": str,       // res.resname (e.g. "LYS", "HEM")
          "het_flag": str,      // res.id[0]: " " or "H_<RESNAME>"
          "resseq": int,        // res.id[1]
          "icode": str,         // res.id[2] (insertion code; usually " ")
          "n_atoms": int,       // len(list(res.get_atoms()))
          "n_unique_atom_ids": int, // len(set(a.id for a in
                                    //   res.get_unpacked_list()))
                                    //   == DockQ list_atoms_per_residue count
                                    //   (altloc-expanded, deduplicated by id)
          "atoms": [
            {
              "name": str,            // atom.id (atom name, e.g. "CA")
              "element": str | null,  // atom.element (e.g. "C", "N")
              "altloc": str,          // atom.get_altloc() (usually " ")
              "coord": [float,float,float],      // float32 widened to float
              "coord_hex": [str,str,str]         // float.hex of each component
            }
          ]
        }
      ]
    }
  ]
}
```

**Why `n_unique_atom_ids` matters:** DockQ's `list_atoms_per_residue` counts
`len(set(a.id for a in res.get_unpacked_list()))` per residue and feeds that
array (`atoms_per_res`) into the `residue_distances` kernel. The atom *coords*
fed to the kernel come from `res.get_atoms()` (count = `n_atoms`). For residues
without alternate locations the two counts are equal; a Rust port must
replicate both semantics (coords from `get_atoms`, per-residue atom *count* from
the deduplicated unpacked-list ids).

---

## (B) Interface dumps — `dumps/interfaces/`

One file per curated case (`<case>.json`). Each replicates `DockQ.main()`'s
chain-mapping search to find the best `chain_map`, then dumps the per-interface
intermediates for that map. Cases and their CLI equivalents:

| case | model | native | flags | best mapping | GlobalDockQ |
|------|-------|--------|-------|--------------|-------------|
| `1A2K`         | `1A2K_r_l_b.model.pdb` | `1A2K_r_l_b.pdb` | — | `BAC:ABC` | 0.653 |
| `1A2K_noalign` | `1A2K_r_l_b.model.pdb` | `1A2K_r_l_b.pdb` | `--no_align` | `BAC:ABC` | 0.653 |
| `dimer_dimer`  | `dimer_dimer.model.pdb` | `dimer_dimer.pdb` | — | `ABLH:ABLH` | 1.000 |
| `model_native` | `model.pdb` | `native.pdb` | `--allowed_mismatches 1` | `AB:AB` | 0.700 |
| `1EXB`         | `1EXB_r_l_b.model.pdb` | `1EXB_r_l_b.pdb` | — | `BACDFHEG:ABDCEGFH` | 0.852 |
| `1EXB_cif`     | `1EXB_r_l_b.model.pdb` | `1EXB.cif.gz` | — | `DH:AE` | 0.775 |
| `6qwn_capri`   | `6qwn-assembly1.cif.gz` | `6qwn-assembly2.cif.gz` | `--capri_peptide` | `AF:BG` | 0.872 |

These match `testdata/*.dockq` exactly. `model_native` uses
`--allowed_mismatches 1` because `model.pdb`/`native.pdb` differ by 1 residue
(per `run_test.sh` line 25). `6qwn_capri` uses `--capri_peptide` (per
`run_test.sh` line 47; golden `testdata/6q2n_peptide.dockq`).

```jsonc
{
  "name": str,                  // case name
  "model": str, "native": str,  // absolute paths
  "kwargs": { ... },            // flags applied (no_align/capri_peptide/
                                //   allowed_mismatches as relevant)
  "error": null | str,          // traceback if the case failed
  "best_chain_map": {nat: mod}, // best native->model chain map (dict)
  "best_mapping_str": str,      // format_mapping_string(best_chain_map),
                                //   "<modelchains>:<nativechains>"
  "best_total_dockq": float,    // sum of per-interface DockQ for best map
  "total_dockq": float,         // == best_total_dockq (from authoritative rerun)
  "n_interfaces": int,          // number of scored native interfaces
  "GlobalDockQ": float,         // total_dockq / n_interfaces
  "interfaces": {
    "<natpair>": {              // e.g. "AB" = native chains A,B (sorted combo)
      "is_het": bool,           // true if either native chain is a het ligand
      "info": { ... },          // full info dict (see below)
      // --- protein-protein interfaces only (is_het == false): ---
      "alignments": [           // one per chain in the pair, native order
        { "seqA": str,          //   model aligned sequence (gaps as "-")
          "matches": str,       //   "|" match, "." mismatch, " " gap
          "seqB": str }         //   native aligned sequence
      ],
      "residue_distances": { ... },        // model (sample) matrix; see below
      "native_residue_distances": { ... }, // native (ref) matrix; same schema
      "fnat_stats": [int,int,int,int],     // [n_shared, n_nonnat, n_native,
                                           //   n_model] from get_fnat_stats
                                           //   on (sample, ref, fnat_threshold)
      "irms": float,            // interface RMSD (Angstrom), == info.iRMSD
      "lrms": float,            // ligand RMSD (Angstrom), == info.LRMSD
      // --- het / small-molecule interfaces only (is_het == true): ---
      "alignments": [ {seqA,matches,seqB}, ... ],  // receptor + ligand align
      "lrms": float | null      // == info.LRMSD (symmetry-corrected)
    }
  }
}
```

### `residue_distances` block

Both `residue_distances` (model/sample) and `native_residue_distances` (ref)
use this schema. The matrix is squared min-atom distances (Angstrom^2), shape
`[n_res_chain1, n_res_chain2]`, exactly as fed to `get_fnat_stats`.

```jsonc
{
  "shape": [int, int],          // [rows, cols]
  "n_lt_25": int,               // count(value < 25.0)  == contacts at < 5 A
  "n_lt_100": int,              // count(value < 100.0) == interface at < 10 A
  "n_lt_4": int,                // count(value < 4.0)   == clashes at < 2 A
  "n_lt_fnat": int,             // count(value < fnat_threshold^2) for this case
  "n_lt_interface": int,        // count(value < interface_threshold^2)
  "fnat_threshold": float,      // 5.0 normally, 4.0 for capri_peptide
  "interface_threshold": float, // 10.0 normally, 8.0 for capri_peptide
  "dtype": str,                 // numpy dtype of the matrix as computed
  "bytes_sha256": str,          // sha256 of matrix.astype(float32).tobytes()
                                //   (row-major / C-order, little-endian)
  "full": [[float,...],...] | null  // nested list iff rows*cols <= 4000,
                                    //   else null (all real cases here: null)
}
```

`n_lt_25 / n_lt_100 / n_lt_4` are fixed 5/10/2 Angstrom references for stable
cross-checking regardless of case. `n_lt_fnat / n_lt_interface` use the case's
actual thresholds (they differ for `--capri_peptide`).

### `info` dict — protein-protein interface (18 keys)

Returned by `run_on_all_native_interfaces` with `low_memory=False`:

| key | type | meaning (units) |
|-----|------|-----------------|
| `DockQ` | float | DockQ score [0,1] |
| `F1` | float | F1 of native contacts |
| `iRMSD` | float | interface RMSD (Angstrom) |
| `LRMSD` | float | ligand RMSD (Angstrom) |
| `fnat` | float | fraction native contacts recovered |
| `fnonnat` | float | fraction model contacts that are non-native |
| `nat_correct` | int | shared (correct) native contacts |
| `nat_total` | int | total native contacts |
| `nonnat_count` | int | non-native model contacts |
| `model_total` | int | total model contacts |
| `clashes` | int | model residue pairs with min-dist < 2 Angstrom |
| `len1` | int | native chain-group-1 residue count |
| `len2` | int | native chain-group-2 residue count |
| `class1` | str | "receptor" or "ligand" (group 1) |
| `class2` | str | "receptor" or "ligand" (group 2) |
| `chain1` | str | model chain id mapped to native pair[0] |
| `chain2` | str | model chain id mapped to native pair[1] |
| `is_het` | bool | always `false` for protein interfaces |
| `chain_map` | dict | the native->model chain map used (diagnostic) |

### `info` dict — het / small-molecule interface (7 keys)

Returned by `calc_sym_corrected_lrmsd` (only the LRMSD-based score):

| key | type | meaning |
|-----|------|---------|
| `DockQ` | float | `dockq_formula(0, 0, LRMSD)` |
| `LRMSD` | float | symmetry-corrected ligand RMSD (Angstrom) |
| `mapping` | dict | best atom isomorphism (sample idx -> ref idx) |
| `is_het` | str | ligand resname (e.g. "HEM") |
| `chain1`, `chain2` | str | model chain ids |
| `chain_map` | dict | native->model chain map (diagnostic) |

(Het interfaces appear only when running with `--small_molecule`; none of the
seven curated cases above are het. The schema is documented for the
`oracle_run.py --small_molecule` differential target.)

---

## (C) CLI differential target — `oracle_run.py`

Reproduces `DockQ.main()`'s chain-mapping search and emits ONE JSON object to
stdout. Use as a stdout-diff target for the Rust CLI.

```
python oracle_run.py <model> <native> [--no_align] [--capri_peptide] \
                     [--mapping M:N] [--allowed_mismatches N] [--small_molecule]
```

Output (stdout, `sort_keys=True`, full precision, trailing newline):

```jsonc
{
  "model": str,
  "native": str,
  "best_mapping_str": str,       // "<modelchains>:<nativechains>"
  "best_dockq": float,           // summed DockQ over best interfaces
  "GlobalDockQ": float,          // best_dockq / number_of_interfaces
  "best_result": { "<natpair>": info_dict, ... }   // info dict as in (B)
}
```

On any exception it prints `{"error": "<traceback>"}` to stdout and exits **1**
(explicit, never silent). Success exits **0**.

The mapping search mirrors `main()` exactly but runs **serially** (n_cpu=1
equivalent) for determinism. Selection is `total_dockq > best_dockq`
(first-wins on ties) over the same order-preserving `chain_maps` iterator
`main()` uses, so the result is identical to `main()`'s parallel reduction. The
`low_memory` (>100 combos) rerun-the-best branch is preserved.

Verified against golden CLI (`testdata/`):
- 1A2K -> `BAC:ABC`, GlobalDockQ 0.6529 (CLI "over 3 native interfaces: 0.653")
- 1EXB -> `BACDFHEG:ABDCEGFH`, GlobalDockQ 0.8516 (16 interfaces, 0.852)
- 6qwn `--capri_peptide` -> `AF:BG`, 0.872; model/native `--allowed_mismatches 1`
  -> `AB:AB`, 0.700; per-field DockQ/iRMSD/LRMSD/fnat/fnonnat/F1/clashes match.

---

## (E) Regenerate everything

Activate the venv in every shell first:

```bash
source /Users/Andre.Teixeira/projects/DockQ/.venv-baseline/bin/activate
```

Then:

```bash
# (A) parse dumps -> dumps/parse/*.json and dumps/parse_het/*.json
python /Users/Andre.Teixeira/projects/dockq-rs/oracle/dump_parse.py

# (B) interface dumps -> dumps/interfaces/*.json
python /Users/Andre.Teixeira/projects/dockq-rs/oracle/dump_interfaces.py

# (C) CLI differential target (examples)
python /Users/Andre.Teixeira/projects/dockq-rs/oracle/oracle_run.py \
    /Users/Andre.Teixeira/projects/DockQ/examples/1A2K_r_l_b.model.pdb \
    /Users/Andre.Teixeira/projects/DockQ/examples/1A2K_r_l_b.pdb

python /Users/Andre.Teixeira/projects/dockq-rs/oracle/oracle_run.py \
    /Users/Andre.Teixeira/projects/DockQ/examples/6qwn-assembly1.cif.gz \
    /Users/Andre.Teixeira/projects/DockQ/examples/6qwn-assembly2.cif.gz --capri_peptide
```

Validate all emitted JSON:

```bash
python - <<'PY'
import json, glob
for f in sorted(glob.glob('/Users/Andre.Teixeira/projects/dockq-rs/oracle/dumps/**/*.json', recursive=True)):
    json.load(open(f))
print("all JSON valid")
PY
```

Scripts are idempotent: rerunning overwrites the dumps in place.
