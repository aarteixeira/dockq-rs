//! End-to-end integration tests on vendored example structures (no external deps —
//! `cargo test` proves the full parse → align → score → mapping-search pipeline anywhere).
//! Expected values are the reference DockQ v2.1.3 results.

use dockq_core::{score_pair, DockQError, RunOptions};

fn data(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

#[test]
fn scores_1a2k_matches_reference() {
    let r = score_pair(
        &data("1A2K_r_l_b.model.pdb"),
        &data("1A2K_r_l_b.pdb"),
        &RunOptions::default(),
    )
    .expect("scoring 1A2K");

    // Optimal model:native chain mapping (tie-break parity with the reference).
    assert_eq!(r.best_mapping_str, "BAC:ABC");
    assert_eq!(r.best_result.len(), 3);
    assert!(
        (r.global_dockq - 0.6529).abs() < 1e-3,
        "GlobalDockQ {} != ~0.653",
        r.global_dockq
    );

    // Per-interface DockQ (keyed by native chain pair) vs the reference.
    let dockq = |k: &str| r.best_result.get(k).map(|x| x.dockq).unwrap_or(f64::NAN);
    assert!((dockq("AB") - 0.994).abs() < 2e-3, "AB {}", dockq("AB"));
    assert!((dockq("AC") - 0.511).abs() < 2e-3, "AC {}", dockq("AC"));
    assert!((dockq("BC") - 0.453).abs() < 2e-3, "BC {}", dockq("BC"));

    // fnat on AB is an exact integer ratio.
    let ab = &r.best_result["AB"];
    assert_eq!(ab.nat_total, 119);
    assert_eq!(ab.nat_correct, 117);
}

#[test]
fn self_comparison_is_perfect() {
    let r = score_pair(
        &data("1A2K_r_l_b.pdb"),
        &data("1A2K_r_l_b.pdb"),
        &RunOptions::default(),
    )
    .expect("self comparison");
    assert!((r.global_dockq - 1.0).abs() < 1e-6, "self GlobalDockQ {}", r.global_dockq);
    assert_eq!(r.best_mapping_str, "ABC:ABC");
}

#[test]
fn mmcif_gz_parses_and_scores() {
    // 1EXB.cif.gz (gzip + mmCIF, chains A/E): self-comparison exercises the full
    // gzip+mmCIF parse → score path and must be a perfect 1.0.
    let r = score_pair(&data("1EXB.cif.gz"), &data("1EXB.cif.gz"), &RunOptions::default())
        .expect("scoring gzipped mmCIF");
    assert!((r.global_dockq - 1.0).abs() < 1e-6, "1EXB self GlobalDockQ {}", r.global_dockq);
}

#[test]
fn small_molecule_flag_errors() {
    let opts = RunOptions {
        small_molecule: true,
        ..RunOptions::default()
    };
    let r = score_pair(
        &data("1A2K_r_l_b.model.pdb"),
        &data("1A2K_r_l_b.pdb"),
        &opts,
    );
    assert!(
        matches!(r, Err(DockQError::SmallMoleculeUnsupported)),
        "small_molecule must hard-error, got {r:?}"
    );
}

#[test]
fn unknown_path_errors_cleanly() {
    let r = score_pair(&data("does_not_exist.pdb"), &data("1A2K_r_l_b.pdb"), &RunOptions::default());
    assert!(matches!(r, Err(DockQError::Io { .. })));
}
