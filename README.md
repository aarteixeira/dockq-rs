# dockq-rs

A Rust reimplementation of the [DockQ](https://github.com/bjornwallner/DockQ) scoring
internals (protein / nucleic-acid core) with a Python wrapper. Same scores as upstream
DockQ v2.1.3, **~20–25× faster** on real workloads, with a parallel batch API and a
Python API that mirrors DockQ's functions (migrate with a one-line import change).

```bash
pip install dockq-rs
```

```text
$ dockq-rs examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --short
Total DockQ over 3 native interfaces: 0.653 with BAC:ABC model:native mapping
DockQ 0.994 iRMSD 0.000 LRMSD 0.000 fnat 0.983 fnonnat 0.008 F1 0.987 clashes 0 mapping BA:AB ...
DockQ 0.511 iRMSD 1.237 LRMSD 6.864 fnat 0.333 fnonnat 0.000 F1 0.500 clashes 0 mapping BC:AC ...
DockQ 0.453 iRMSD 2.104 LRMSD 8.131 fnat 0.500 fnonnat 0.107 F1 0.641 clashes 0 mapping AC:BC ...
```

## Why

DockQ's chain-mapping search over large homomers is slow (the 8-chain 1EXB example explores
576 mappings × 16 interfaces). dockq-rs moves the whole pipeline — PDB/mmCIF parsing,
sequence alignment, the residue-distance kernel, Kabsch superposition, the mapping search,
and batch scoring — into Rust, parallelised with Rayon, while reproducing upstream's results
to f32 precision.

## Status & scope

- **Implemented:** protein and nucleic-acid scoring; all mapping modes (`--mapping`,
  wildcards, `:NATIVE`), `--capri_peptide`, `--no_align`, `--allowed_mismatches`; single +
  batch; PDB and mmCIF, gzip-aware.
- **Deferred:** small-molecule (HETATM) symmetry-corrected LRMSD. `--small_molecule`
  **raises a clear error** rather than silently falling back — see *No silent failures*.
- **Not implemented:** `--optDockQF1` (errors explicitly).

## Performance

Apple Silicon, 12 cores, vs DockQ v2.1.3. Full numbers + caveats in [`bench/RESULTS.md`](bench/RESULTS.md).

| workload | reference | dockq-rs | speedup |
|---|---|---|---|
| 1EXB 576-combo search (16 interfaces) | 4.94 s | 0.21 s | **23.6×** |
| batch: 64 models vs 1 native (per pair) | 64.5 ms | 3.0 ms | **21.3×** |
| 1A2K dimer (full CLI) | 1.92 s | 0.025 s | 77.7× ¹ |

¹ Inflated by Python `multiprocessing` pool startup on a tiny input; the compute-bound
figures (23.6× / 21.3×) are the honest headline.

## Correctness

Validated against the Python/Cython reference as the oracle (see [`tests/`](tests)):

- **Differential, 13 cases** (dimers, 1EXB multimer, mmCIF+gzip, `--no_align`,
  `--allowed_mismatches`, `--capri_peptide`, self-comparisons, 4 mapping modes): worst
  deviation on any score **4.84e-05**; integer contact counts **exact**; mapping / class /
  chain assignment **exact**.
- **Golden parity:** reproduces upstream `testdata/*.dockq` byte-for-byte (8/8).
- **API parity:** `load_PDB` / `run_on_all_native_interfaces` return the same values as the
  reference (verified by running the same logic through both and comparing).
- **30 tests** in `dockq-core` (25 unit + 5 end-to-end integration on vendored data;
  geometry bit-exact to the Cython kernel; alignment fuzzed on 35k+ pairs against Biopython
  with 0 mismatches; parser exact to the f32 bit pattern on all example files).

## Install

```bash
pip install dockq-rs                 # prebuilt wheels for Linux/macOS/Windows
```

From source (needs a Rust toolchain + maturin):

```bash
cd dockq-py && maturin develop --release
```

This provides the `dockq_rs` Python module. The standalone CLI binary is
`cargo build --release -p dockq-cli` (`target/release/dockq-rs`).

## Usage

### CLI (drop-in flags)

```bash
dockq-rs model.pdb native.pdb                       # long output
dockq-rs model.pdb native.pdb --short
dockq-rs model.cif.gz native.cif.gz --capri_peptide
dockq-rs model.pdb native.pdb --mapping 'ABC:BCA'   # or ':BC', 'A*:W*', --no_align, --allowed_mismatches N
dockq-rs model.pdb native.pdb --json out.json
```

### CLI — batch (parallel in Rust)

```bash
# rank many models against one native (TSV, sorted by GlobalDockQ)
dockq-rs batch --native native.pdb --models m1.pdb m2.pdb m3.pdb --sort
dockq-rs batch --native native.pdb --models_dir ./models --sort        # every .pdb/.cif in a dir
dockq-rs batch --native native.pdb --models_from list.txt              # paths, one per line

# arbitrary (model, native) pairs, "model native" per line
dockq-rs batch --pairs_from pairs.txt --format json -o results.json
```

Scoring flags (`--capri_peptide`, `--no_align`, `--allowed_mismatches`, `--mapping`,
`--n_cpu`) apply to every job. Output is TSV (default) or `--format json`. A failed job is
reported as an explicit `error:` row (never silently dropped) and the run exits non-zero if
any job failed.

### Python — migrating from DockQ (one-line import change)

The function signatures and return shapes match the reference; change
`from DockQ.DockQ import ...` to `from dockq_rs import ...`:

```python
from dockq_rs import load_PDB, run_on_all_native_interfaces

model  = load_PDB("model.pdb")
native = load_PDB("native.pdb")
result, total_dockq = run_on_all_native_interfaces(
    model, native, chain_map={"A": "A", "B": "B"}   # native -> model
)
# result[("A","B")]["DockQ"], ["iRMSD"], ["LRMSD"], ["fnat"], ... (same keys as upstream)
```

### Python — new ergonomic + batch API (parallel in Rust)

```python
import dockq_rs

# single
r = dockq_rs.score("model.pdb", "native.pdb")          # dict: GlobalDockQ, best_mapping_str, best_result, ...

# one native, many models (model ranking) — parallel
outcomes = dockq_rs.score_one_vs_many("native.pdb", ["m1.pdb", "m2.pdb", ...])

# arbitrary (model, native) pairs — parallel
outcomes = dockq_rs.score_pairs([("m1.pdb", "n1.pdb"), ("m2.pdb", "n2.pdb")])
# each outcome: {"model","native","ok",  "result"|"error"}   (errors surfaced per job)
```

## No silent failures

Per design, every failure mode is explicit:

- Structure format is detected by **content** (mmCIF `_atom_site` loop vs PDB ATOM records);
  on failure it errors — no "try PDB, except try mmCIF" fallback.
- There is **no** Cython→Python style compute fallback.
- `--small_molecule` and `--optDockQF1` raise typed errors instead of silently doing
  something else.
- Batch jobs report per-job errors in the outcome rather than dropping failures.

## Architecture

```
dockq-core/   Rust library: model, parser (PDB/mmCIF), align (NW), geometry
              (distances/fnat/Kabsch/RMSD), dockq (calc_DockQ), mapping (search), batch.
dockq-cli/    `dockq-rs` binary (output matches upstream, header aside).
dockq-py/     PyO3 extension `dockq_rs` + `dockq_rs`/`DockQ` Python packages (maturin).
oracle/       Python harness dumping reference intermediates/results as JSON ground truth.
tests/        differential.py, golden.sh, dropin_compat.py.
bench/        single.py, batch_*.py, RESULTS.md.
```

Coordinates and geometry are f32 to match Biopython's storage and the Cython kernel.

## Running tests

`cargo test` is self-contained (it uses small vendored example structures, so it passes
on any machine with no external setup):

```bash
cargo test -p dockq-core        # 25 unit + 5 integration tests, no external deps
```

The differential / golden / drop-in tests compare against a reference DockQ checkout.
They auto-discover a sibling `../DockQ`, or point them anywhere with env vars (nothing is
hardcoded):

```bash
export DOCKQ_REPO=/path/to/DockQ                # reference checkout (examples/ + testdata/)
export DOCKQ_PYTHON=$DOCKQ_REPO/.venv/bin/python # python with reference DockQ installed
python3 tests/differential.py                    # vs Python oracle (13 cases)
bash    tests/golden.sh                           # byte-parity vs testdata/*.dockq
```

## Credit

Reimplements the algorithm of DockQ (Mirabello & Wallner; Basu & Wallner). This is an
independent performance-oriented port of the protein/NA scoring path.
