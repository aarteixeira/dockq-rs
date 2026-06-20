"""Drop-in compatibility shim for the original `DockQ` package, backed by the Rust
`dockq_rs` extension. Existing code that does

    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

keeps working unchanged.
"""
from . import DockQ  # noqa: F401
