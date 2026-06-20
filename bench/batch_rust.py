"""Rust batch: score N models against one native via dockq_rs.score_one_vs_many (one
process, parallel in Rust). Same N as bench/batch_ref.py for a like-for-like comparison."""
import sys
import time

import dockq_rs


def main():
    n = int(sys.argv[1])
    model = sys.argv[2]
    native = sys.argv[3]
    models = [model] * n
    t0 = time.perf_counter()
    res = dockq_rs.score_one_vs_many(native, models)
    dt = time.perf_counter() - t0
    assert all(x["ok"] for x in res), "a job failed"
    print(f"{dt:.4f}")


if __name__ == "__main__":
    main()
