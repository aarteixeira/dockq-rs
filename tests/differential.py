#!/usr/bin/env python3
"""Differential correctness harness: Rust dockq-rs vs the Python/Cython reference oracle.

For every case we run BOTH:
  - the reference  : oracle/oracle_run.py (real DockQ v2.1.3 in .venv-baseline)
  - the Rust port  : target/release/dockq-rs --diff-json
and compare every reported quantity. Integer contact counts must match EXACTLY (the
distance kernel is bit-exact); RMSDs/DockQ are allowed a tight f32/SVD tolerance; the
chosen mapping, chain assignment and class labels must match exactly.

No silent failures: any mismatch is printed in full and the script exits non-zero.
"""
import json
import subprocess
import sys
from pathlib import Path

DOCKQ = Path("/Users/Andre.Teixeira/projects/DockQ")
EX = DOCKQ / "examples"
PYBIN = DOCKQ / ".venv-baseline/bin/python"
ORACLE = Path("/Users/Andre.Teixeira/projects/dockq-rs/oracle/oracle_run.py")
RUST = Path("/Users/Andre.Teixeira/projects/dockq-rs/target/release/dockq-rs")

# Tolerances
RMSD_ABS = 1.0e-3   # iRMSD / LRMSD / DockQ: f32 geometry + different SVD backend
RATIO_ABS = 1.0e-6  # fnat / fnonnat / F1: derived from exact integer counts
EXACT_INT = ["nat_correct", "nat_total", "nonnat_count", "model_total", "clashes", "len1", "len2"]
EXACT_STR = ["class1", "class2", "chain1", "chain2"]
FLOAT_RMSD = ["DockQ", "iRMSD", "LRMSD"]
FLOAT_RATIO = ["fnat", "fnonnat", "F1"]

# (name, model, native, extra-flags)
CASES = [
    ("1A2K", EX / "1A2K_r_l_b.model.pdb", EX / "1A2K_r_l_b.pdb", []),
    ("1A2K_noalign", EX / "1A2K_r_l_b.model.pdb", EX / "1A2K_r_l_b.pdb", ["--no_align"]),
    ("1A2K_self", EX / "1A2K_r_l_b.pdb", EX / "1A2K_r_l_b.pdb", []),
    ("dimer_dimer", EX / "dimer_dimer.model.pdb", EX / "dimer_dimer.pdb", []),
    ("dimer_dimer_self", EX / "dimer_dimer.pdb", EX / "dimer_dimer.pdb", []),
    ("model_native", EX / "model.pdb", EX / "native.pdb", ["--allowed_mismatches", "1"]),
    ("1EXB", EX / "1EXB_r_l_b.model.pdb", EX / "1EXB_r_l_b.pdb", []),
    ("1EXB_AB.BA", EX / "1EXB_r_l_b.model.pdb", EX / "1EXB_r_l_b.pdb", ["--mapping", "AB*:BA*"]),
    ("1EXB_.ABC", EX / "1EXB_r_l_b.model.pdb", EX / "1EXB_r_l_b.pdb", ["--mapping", ":ABC"]),
    ("1EXB_full", EX / "1EXB_r_l_b.model.pdb", EX / "1EXB_r_l_b.pdb", ["--mapping", "ABCDEFGH:BADCFEHG"]),
    ("1EXB_cif", EX / "1EXB_r_l_b.model.pdb", EX / "1EXB.cif.gz", ["--mapping", "DH:AE"]),
    ("1EXB_self", EX / "1EXB_r_l_b.pdb", EX / "1EXB_r_l_b.pdb", []),
    ("6qwn_capri", EX / "6qwn-assembly1.cif.gz", EX / "6qwn-assembly2.cif.gz", ["--capri_peptide"]),
]


def run_json(cmd):
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if p.returncode != 0:
        return None, f"exit {p.returncode}: {p.stderr.strip()[:500]}"
    try:
        return json.loads(p.stdout), None
    except json.JSONDecodeError as e:
        return None, f"bad JSON: {e}; stdout head: {p.stdout[:300]}"


def run_oracle(model, native, flags):
    return run_json([PYBIN, ORACLE, model, native, *flags])


def run_rust(model, native, flags):
    return run_json([RUST, model, native, "--diff_json", *flags])


def cmp_interface(name, pair, o, r):
    errs = []
    for k in EXACT_INT:
        if int(o[k]) != int(r[k]):
            errs.append(f"  [{pair}] {k}: oracle={o[k]} rust={r[k]} (must be exact)")
    for k in EXACT_STR:
        if str(o[k]) != str(r[k]):
            errs.append(f"  [{pair}] {k}: oracle={o[k]!r} rust={r[k]!r}")
    for k in FLOAT_RATIO:
        if abs(float(o[k]) - float(r[k])) > RATIO_ABS:
            errs.append(f"  [{pair}] {k}: oracle={o[k]:.9f} rust={r[k]:.9f} d={abs(o[k]-r[k]):.2e}")
    for k in FLOAT_RMSD:
        if abs(float(o[k]) - float(r[k])) > RMSD_ABS:
            errs.append(f"  [{pair}] {k}: oracle={o[k]:.6f} rust={r[k]:.6f} d={abs(o[k]-r[k]):.2e}")
    return errs


def cmp_case(name, o, r):
    errs = []
    # Chosen mapping must match exactly (tie-break parity).
    if o.get("best_mapping_str") != r.get("best_mapping_str"):
        errs.append(f"  best_mapping_str: oracle={o.get('best_mapping_str')!r} rust={r.get('best_mapping_str')!r}")
    if abs(float(o["GlobalDockQ"]) - float(r["GlobalDockQ"])) > RMSD_ABS:
        errs.append(f"  GlobalDockQ: oracle={o['GlobalDockQ']:.6f} rust={r['GlobalDockQ']:.6f}")
    ob, rb = o["best_result"], r["best_result"]
    if set(ob.keys()) != set(rb.keys()):
        errs.append(f"  interface keys differ: oracle={sorted(ob)} rust={sorted(rb)}")
        return errs
    for pair in ob:
        errs += cmp_interface(name, pair, ob[pair], rb[pair])
    return errs


def main():
    if not RUST.exists():
        print(f"FATAL: Rust binary not built at {RUST}. Run: cargo build --release -p dockq-cli")
        sys.exit(2)
    total = 0
    failed = 0
    max_rmsd_diff = 0.0
    for name, model, native, flags in CASES:
        total += 1
        o, oerr = run_oracle(model, native, flags)
        r, rerr = run_rust(model, native, flags)
        if oerr:
            print(f"FAIL {name}: oracle error: {oerr}")
            failed += 1
            continue
        if rerr:
            print(f"FAIL {name}: rust error: {rerr}")
            failed += 1
            continue
        if "error" in o:
            print(f"FAIL {name}: oracle returned error: {o['error'][:300]}")
            failed += 1
            continue
        errs = cmp_case(name, o, r)
        # track worst rmsd diff for reporting
        for pair in o["best_result"]:
            if pair in r["best_result"]:
                for k in FLOAT_RMSD:
                    max_rmsd_diff = max(max_rmsd_diff, abs(float(o["best_result"][pair][k]) - float(r["best_result"][pair][k])))
        if errs:
            failed += 1
            print(f"FAIL {name} ({len(errs)} diffs):")
            for e in errs:
                print(e)
        else:
            nif = len(o["best_result"])
            print(f"PASS {name}: GlobalDockQ={o['GlobalDockQ']:.3f} mapping={o['best_mapping_str']} ({nif} interfaces)")
    print(f"\n{total - failed}/{total} cases passed. Worst RMSD/DockQ abs diff vs oracle: {max_rmsd_diff:.2e}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
