"""Single-pair latency: reference DockQ CLI vs Rust dockq-rs CLI (median of K runs,
full process incl. parse + chain-mapping search). The 1EXB search case (576 chain-mapping
combinations, 16 interfaces) is the "calculations take too long" headline."""
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

# Portable discovery (override with DOCKQ_REPO / DOCKQ_PYTHON / DOCKQ_RS_BIN).
REPO = Path(__file__).resolve().parent.parent
RUST = str(os.environ.get("DOCKQ_RS_BIN", REPO / "target" / "release" / "dockq-rs"))
DOCKQ = os.environ.get("DOCKQ_REPO") or str(REPO.parent / "DockQ")
if not (Path(DOCKQ) / "examples").is_dir():
    sys.exit("Set DOCKQ_REPO to your reference DockQ checkout (with examples/).")
PYBIN = os.environ.get("DOCKQ_PYTHON", f"{DOCKQ}/.venv-baseline/bin/python")

K = int(sys.argv[1]) if len(sys.argv) > 1 else 7

CASES = [
    ("1A2K (dimer)", ["examples/1A2K_r_l_b.model.pdb", "examples/1A2K_r_l_b.pdb", "--short"]),
    ("dimer_dimer", ["examples/dimer_dimer.model.pdb", "examples/dimer_dimer.pdb", "--short"]),
    ("1EXB (576-combo search)", ["examples/1EXB_r_l_b.model.pdb", "examples/1EXB_r_l_b.pdb", "--short"]),
    ("1EXB (fixed mapping)", ["examples/1EXB_r_l_b.model.pdb", "examples/1EXB_r_l_b.pdb", "--short", "--mapping", "ABCDEFGH:BADCFEHG"]),
]


def timeit(cmd):
    ts = []
    for _ in range(K):
        t0 = time.perf_counter()
        p = subprocess.run(cmd, cwd=DOCKQ, capture_output=True)
        ts.append(time.perf_counter() - t0)
        if p.returncode != 0:
            return None
    return statistics.median(ts)


print(f"single-pair latency (median of {K}, full CLI incl. startup):")
print(f"{'case':28s} {'python':>11s} {'rust':>10s} {'speedup':>9s}")
for name, args in CASES:
    pt = timeit([PYBIN, "-m", "DockQ.DockQ"] + args)
    rt = timeit([RUST] + args)
    if pt is None or rt is None:
        print(f"{name:28s}  ERROR")
        continue
    print(f"{name:28s} {pt*1000:9.1f}ms {rt*1000:8.1f}ms {pt/rt:8.1f}x")
