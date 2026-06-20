"""dockq_rs — Rust reimplementation of the DockQ scoring internals (protein/NA core).

Public API:
  * load_PDB, run_on_all_native_interfaces  — drop-in compatible with reference DockQ.
  * score, score_one_vs_many, score_pairs    — new ergonomic + batch API (parallel in Rust).
"""
from ._dockq_rs import (  # noqa: F401
    Structure,
    load_PDB,
    run_on_all_native_interfaces,
    score,
    score_one_vs_many,
    score_pairs,
    version,
)

__all__ = [
    "Structure",
    "load_PDB",
    "run_on_all_native_interfaces",
    "score",
    "score_one_vs_many",
    "score_pairs",
    "version",
]
__version__ = version()
