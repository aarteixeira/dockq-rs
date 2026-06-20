//! Chain grouping + mapping enumeration + parallel search (task #6, integration owner).
//! Ports `group_chains`, `product_without_dupl` (exact enumeration order), `get_all_chain_maps`,
//! `run_on_all_native_interfaces`, and the multimer search. Tie-break on equal total DockQ
//! is "first in enumeration order wins" (strict `>`), preserved under Rayon via a
//! deterministic argmax over (dockq, enumeration index).
