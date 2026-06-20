#!/usr/bin/env python
"""Ground-truth oracle: CLI differential target.

Reproduces DockQ.main()'s chain-mapping search and prints ONE JSON object with
the best mapping and full per-interface results. Intended as a stdout-diff
target for a Rust port.

Usage:
  python oracle_run.py <model> <native> [--no_align] [--capri_peptide]
                       [--mapping M:N] [--allowed_mismatches N] [--small_molecule]

Output (stdout, sorted keys, full numeric precision):
  {"model":..., "native":..., "best_mapping_str":..., "best_dockq":...,
   "GlobalDockQ":..., "best_result": {"<natpair>": info_dict, ...}}

On any exception: prints {"error": "<traceback>"} to stdout and exits 1.

The mapping search mirrors main() exactly but runs serially (n_cpu=1 equivalent)
for determinism. Selection is total_dockq > best_dockq (first-wins ties), which
is identical to main()'s reduction over the order-preserving chain_maps iterator.
"""
import argparse
import itertools
import json
import sys
import traceback

import numpy as np

from DockQ.DockQ import (
    count_chain_combinations,
    format_mapping,
    format_mapping_string,
    get_all_chain_maps,
    group_chains,
    load_PDB,
    run_on_all_native_interfaces,
)


def parse_args():
    p = argparse.ArgumentParser(description="DockQ oracle CLI differential target")
    p.add_argument("model")
    p.add_argument("native")
    p.add_argument("--no_align", action="store_true")
    p.add_argument("--capri_peptide", action="store_true")
    p.add_argument("--mapping", default=None)
    p.add_argument("--allowed_mismatches", type=int, default=0)
    p.add_argument("--small_molecule", action="store_true")
    return p.parse_args()


def json_safe(obj):
    """Recursively cast numpy scalars/arrays to native types for JSON."""
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    return obj


def run(args):
    initial_mapping, model_chains, native_chains = format_mapping(
        args.mapping, args.small_molecule
    )
    model_structure = load_PDB(
        args.model, chains=model_chains or [], small_molecule=args.small_molecule
    )
    native_structure = load_PDB(
        args.native, chains=native_chains or [], small_molecule=args.small_molecule
    )
    model_chains = (
        [c.id for c in model_structure] if not model_chains else model_chains
    )
    native_chains = (
        [c.id for c in native_structure] if not native_chains else native_chains
    )
    if len(model_chains) < 2 or len(native_chains) < 2:
        raise ValueError("Need at least two chains in the two inputs")

    model_chains_to_combo = [
        mc for mc in model_chains if mc not in initial_mapping.values()
    ]
    native_chains_to_combo = [
        nc for nc in native_chains if nc not in initial_mapping.keys()
    ]
    chain_clusters, reverse_map = group_chains(
        model_structure,
        native_structure,
        model_chains_to_combo,
        native_chains_to_combo,
        args.allowed_mismatches,
    )
    chain_maps = get_all_chain_maps(
        chain_clusters,
        initial_mapping,
        reverse_map,
        model_chains_to_combo,
        native_chains_to_combo,
    )
    num_chain_combinations = count_chain_combinations(chain_clusters)
    chain_maps, chain_maps_ = itertools.tee(chain_maps)
    low_memory = num_chain_combinations > 100

    best_dockq = -1
    best_result = None
    best_mapping = None

    if num_chain_combinations > 1:
        for chain_map in chain_maps_:
            result_this_mapping, total_dockq = run_on_all_native_interfaces(
                model_structure,
                native_structure,
                chain_map=chain_map,
                no_align=args.no_align,
                capri_peptide=args.capri_peptide,
                low_memory=low_memory,
            )
            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_result = result_this_mapping
                best_mapping = chain_map
        if low_memory:
            best_result, best_dockq = run_on_all_native_interfaces(
                model_structure,
                native_structure,
                chain_map=best_mapping,
                no_align=args.no_align,
                capri_peptide=args.capri_peptide,
                low_memory=False,
            )
    else:
        best_mapping = next(chain_maps)
        best_result, best_dockq = run_on_all_native_interfaces(
            model_structure,
            native_structure,
            chain_map=best_mapping,
            no_align=args.no_align,
            capri_peptide=args.capri_peptide,
            low_memory=False,
        )

    if not best_result:
        raise RuntimeError(
            "Could not find interfaces in the native model. "
            "Check inputs or select chains with --mapping."
        )

    out = {
        "model": args.model,
        "native": args.native,
        "best_mapping_str": format_mapping_string(best_mapping),
        "best_dockq": float(best_dockq),
        "GlobalDockQ": float(best_dockq) / len(best_result),
        "best_result": {k: json_safe(v) for k, v in best_result.items()},
    }
    return out


def main():
    args = parse_args()
    try:
        out = run(args)
    except Exception:
        json.dump({"error": traceback.format_exc()}, sys.stdout, sort_keys=True)
        sys.stdout.write("\n")
        sys.exit(1)
    json.dump(out, sys.stdout, sort_keys=True)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
