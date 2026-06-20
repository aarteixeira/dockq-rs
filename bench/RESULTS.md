# Benchmark results

Hardware: Apple Silicon, 12 cores. Rust built `--release`. Reference: DockQ v2.1.3
(Python 3.12 + Cython) in `.venv-baseline`. Reproduce with `bench/single.py` and
`bench/batch_*.py` (see below).

## Speed

### Single pair — full CLI, median of 7 runs (parse + chain-mapping search + output)

| case | reference (python) | dockq-rs | speedup |
|---|---|---|---|
| 1A2K (dimer) | 1923.8 ms | 24.8 ms | **77.7×** |
| dimer_dimer (4 interfaces) | 685.6 ms | 35.8 ms | **19.2×** |
| 1EXB (576-combo search, 16 interfaces) | 4942.9 ms | 209.3 ms | **23.6×** |
| 1EXB (fixed mapping) | 1064.0 ms | 57.6 ms | **18.5×** |

### Batch — 64 models vs 1 native (in-process; no per-call interpreter startup)

| | total | per pair | speedup |
|---|---|---|---|
| reference (sequential, in-process) | 4.128 s | 64.5 ms | — |
| dockq-rs `score_one_vs_many` (parallel) | 0.194 s | 3.0 ms | **21.3×** |

## Honest caveats

- The single-pair numbers are **full process** times (the way a user actually runs the
  tool), so the Python column includes interpreter startup **and** its `multiprocessing`
  pool spin-up (`parallelbar.progress_map`). For tiny inputs that overhead dominates — the
  1A2K **77.7×** is mostly Python paying pool-setup cost for a 3-combination search, not a
  pure-compute gap. The most representative *compute-bound* figure is **1EXB search:
  4.9 s → 0.21 s (23.6×)** and the **in-process batch: 64.5 → 3.0 ms/pair (21.3×)**, which
  exclude per-call startup on both sides.
- The batch test reuses one model path ×64 as a stand-in for 64 distinct models; per-call
  work (parse + search) is identical to distinct files, but OS file-cache effects are
  removed from the comparison (favouring neither side equally — both re-parse each call).
- `residue_distances` micro-benchmark (the dominant kernel, in isolation): 300×300
  residues in **0.36 ms** parallel (5.8× over serial via Rayon).

## Accuracy

Differential vs the Python/Cython oracle across 13 cases (dimers, the 16-interface 1EXB
multimer, mmCIF+gzip, `--no_align`, `--allowed_mismatches`, `--capri_peptide`,
self-comparisons, and four mapping modes):

- **Worst absolute deviation on any DockQ / iRMSD / LRMSD value: 4.84e-05** (f32 geometry
  + a different SVD backend than NumPy's LAPACK).
- **Integer contact counts** (nat_correct, nat_total, nonnat_count, model_total, clashes,
  len1, len2) match **exactly** — the distance kernel is bit-identical to the Cython one.
- Chosen chain mapping, receptor/ligand class, and chain assignment match **exactly**
  (tie-break parity).
- Golden `testdata/*.dockq` CLI output reproduced **byte-for-byte** (header aside): 8/8.
- Drop-in API parity: identical user code in both venvs yields identical results
  (1A2K `{"A":"A","B":"B"}` → 0.9425, matching the reference exactly).

## Reproduce

```bash
# single-pair latency (K=7)
python3 bench/single.py 7

# batch (N=64)
EX=/path/to/DockQ/examples
.venv-baseline/bin/python bench/batch_ref.py 64 $EX/1A2K_r_l_b.model.pdb $EX/1A2K_r_l_b.pdb
dockq-rs-venv/bin/python  bench/batch_rust.py 64 $EX/1A2K_r_l_b.model.pdb $EX/1A2K_r_l_b.pdb
```
