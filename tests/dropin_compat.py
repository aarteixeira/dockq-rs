"""Drop-in API parity: identical user code, run in BOTH the reference venv and the Rust
`dockq_rs` venv, comparing `load_PDB` + `run_on_all_native_interfaces` outputs.

Run in each venv with a tag, then diff the two JSON files:
    .venv-baseline/bin/python tests/dropin_compat.py ref  > /tmp/dropin_ref.json
    dockq-rs-venv/bin/python  tests/dropin_compat.py rust > /tmp/dropin_rust.json
    python tests/dropin_compat.py diff /tmp/dropin_ref.json /tmp/dropin_rust.json
"""
import json
import os
import sys
from pathlib import Path


def _examples_dir():
    env = os.environ.get("DOCKQ_REPO")
    if env:
        return f"{env}/examples"
    sibling = Path(__file__).resolve().parent.parent.parent / "DockQ" / "examples"
    if sibling.is_dir():
        return str(sibling)
    sys.exit("Set DOCKQ_REPO to your reference DockQ checkout (with examples/).")


EX = _examples_dir()

# (name, model, native, chain_map native->model)
CASES = [
    ("1A2K_AB", f"{EX}/1A2K_r_l_b.model.pdb", f"{EX}/1A2K_r_l_b.pdb", {"A": "A", "B": "B"}),
    ("1A2K_BA", f"{EX}/1A2K_r_l_b.model.pdb", f"{EX}/1A2K_r_l_b.pdb", {"A": "B", "B": "A"}),
    ("1A2K_ABC", f"{EX}/1A2K_r_l_b.model.pdb", f"{EX}/1A2K_r_l_b.pdb", {"A": "B", "B": "A", "C": "C"}),
    ("dimer", f"{EX}/dimer_dimer.model.pdb", f"{EX}/dimer_dimer.pdb", {"A": "A", "B": "B", "L": "L", "H": "H"}),
]

COMMON = ["DockQ", "F1", "iRMSD", "LRMSD", "fnat", "fnonnat", "nat_correct", "nat_total",
          "nonnat_count", "model_total", "clashes", "len1", "len2", "class1", "class2",
          "chain1", "chain2"]


def collect():
    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
    out = {}
    for name, m, n, cmap in CASES:
        model = load_PDB(m)
        native = load_PDB(n)
        result, total = run_on_all_native_interfaces(model, native, chain_map=cmap)
        ifaces = {}
        for k, v in result.items():
            key = f"{k[0]}{k[1]}"
            ifaces[key] = {kk: v[kk] for kk in COMMON if kk in v}
        out[name] = {"total": total, "interfaces": ifaces}
    print(json.dumps(out, sort_keys=True))


def diff(ref_path, rust_path):
    ref = json.load(open(ref_path))
    rust = json.load(open(rust_path))
    fails = 0
    for name in ref:
        if name not in rust:
            print(f"FAIL {name}: missing in rust")
            fails += 1
            continue
        if abs(ref[name]["total"] - rust[name]["total"]) > 1e-3:
            print(f"FAIL {name}: total {ref[name]['total']} vs {rust[name]['total']}")
            fails += 1
        ri, ti = ref[name]["interfaces"], rust[name]["interfaces"]
        if set(ri) != set(ti):
            print(f"FAIL {name}: interface keys {sorted(ri)} vs {sorted(ti)}")
            fails += 1
            continue
        for pair in ri:
            for k in COMMON:
                rv, tv = ri[pair].get(k), ti[pair].get(k)
                if isinstance(rv, str) or isinstance(tv, str):
                    if rv != tv:
                        print(f"FAIL {name}/{pair}/{k}: {rv!r} vs {tv!r}")
                        fails += 1
                elif rv is None or tv is None:
                    if rv != tv:
                        print(f"FAIL {name}/{pair}/{k}: {rv} vs {tv}")
                        fails += 1
                else:
                    tol = 1e-6 if k in ("fnat", "fnonnat", "F1") else 1e-3
                    if abs(float(rv) - float(tv)) > tol:
                        print(f"FAIL {name}/{pair}/{k}: {rv} vs {tv} (d={abs(rv-tv):.2e})")
                        fails += 1
        if not any(False for _ in []):
            print(f"PASS {name}: total={ref[name]['total']:.4f} ({len(ri)} interfaces)")
    print(f"\ndrop-in parity: {'OK' if fails == 0 else f'{fails} FAILURES'}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    if sys.argv[1] == "diff":
        diff(sys.argv[2], sys.argv[3])
    else:
        collect()
