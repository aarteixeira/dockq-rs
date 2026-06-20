"""`DockQ.DockQ` drop-in surface, re-exported from the Rust `dockq_rs` extension.

The reference's documented public API is preserved with identical names, arguments and
return shapes:

    model  = load_PDB("model.pdb")
    native = load_PDB("native.pdb")
    result, total_dockq = run_on_all_native_interfaces(
        model, native, chain_map={"A": "A", "B": "B"}   # native -> model
    )

`result` is a dict keyed by (native_chain1, native_chain2) tuples, each value a dict with
the same fields as upstream (DockQ, F1, iRMSD, LRMSD, fnat, nat_correct, nat_total,
fnonnat, nonnat_count, model_total, clashes, len1, len2, class1, class2, chain1, chain2,
is_het). Small-molecule scoring is not implemented in this build and raises (no silent
fallback).

New helpers (not in upstream): `score`, `score_one_vs_many`, `score_pairs`.
"""
from dockq_rs import (  # noqa: F401
    load_PDB,
    run_on_all_native_interfaces,
    score,
    score_one_vs_many,
    score_pairs,
    version,
)

__all__ = [
    "load_PDB",
    "run_on_all_native_interfaces",
    "score",
    "score_one_vs_many",
    "score_pairs",
    "version",
]
