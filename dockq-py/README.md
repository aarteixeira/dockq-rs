# dockq-rs

Fast Rust reimplementation of the [DockQ](https://github.com/bjornwallner/DockQ)
protein / nucleic-acid docking scorer, with a Python API. Same scores as DockQ v2.1.3
(validated to f32 precision against the reference), **~20–25× faster**, with parallel
batch scoring.

```bash
pip install dockq-rs
```

## Usage

```python
import dockq_rs

# single model vs native
r = dockq_rs.score("model.pdb", "native.pdb")
print(r["GlobalDockQ"], r["best_mapping_str"])

# rank many models against one native (parallel in Rust)
outcomes = dockq_rs.score_one_vs_many("native.pdb", ["m1.pdb", "m2.pdb", "m3.pdb"])

# arbitrary (model, native) pairs (parallel)
outcomes = dockq_rs.score_pairs([("m1.pdb", "n1.pdb"), ("m2.pdb", "n2.pdb")])
```

## Migrating from DockQ

The function-level API mirrors the reference. Change:

```python
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
```

to:

```python
from dockq_rs import load_PDB, run_on_all_native_interfaces
```

`load_PDB(path, chains=[], small_molecule=False, n_model=0)` and
`run_on_all_native_interfaces(model, native, chain_map=..., no_align=False,
capri_peptide=False)` keep the same arguments and return shapes — the per-interface dict
has the same keys (`DockQ`, `F1`, `iRMSD`, `LRMSD`, `fnat`, `fnonnat`, `nat_correct`,
`nat_total`, `nonnat_count`, `model_total`, `clashes`, `len1`, `len2`, `class1`, `class2`,
`chain1`, `chain2`).

## Scope

Protein and nucleic-acid scoring, all mapping modes, `capri_peptide`, `no_align`,
`allowed_mismatches`. Small-molecule (HETATM) symmetry-corrected scoring is **not**
implemented in this build and raises a clear error (no silent fallback).

## Links

Source, standalone CLI, benchmarks, and correctness details:
<https://github.com/aarteixeira/dockq-rs>. MIT licensed.
