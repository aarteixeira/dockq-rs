"""Reference (Python) batch baseline: score N model/native pairs in-process, sequentially
(the reference has no parallel batch API; this is the realistic way to batch it, without
paying per-call interpreter startup). Replicates main()'s chain-mapping search per pair."""
import sys
import time

from DockQ.DockQ import (
    get_all_chain_maps,
    group_chains,
    load_PDB,
    run_on_all_native_interfaces,
)


def score_one(model_path, native_path):
    model = load_PDB(model_path)
    native = load_PDB(native_path)
    model_chains = [c.id for c in model]
    native_chains = [c.id for c in native]
    chain_clusters, reverse_map = group_chains(model, native, model_chains, native_chains, 0)
    chain_maps = get_all_chain_maps(chain_clusters, {}, reverse_map, model_chains, native_chains)
    best = -1.0
    for cm in chain_maps:
        _res, total = run_on_all_native_interfaces(model, native, chain_map=cm)
        if total > best:
            best = total
    return best


def main():
    n = int(sys.argv[1])
    model = sys.argv[2]
    native = sys.argv[3]
    t0 = time.perf_counter()
    for _ in range(n):
        score_one(model, native)
    dt = time.perf_counter() - t0
    print(f"{dt:.4f}")


if __name__ == "__main__":
    main()
