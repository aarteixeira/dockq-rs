//! dockq-rs CLI. Drop-in-compatible single-pair scoring (output matches the reference
//! DockQ, header aside), plus `--diff-json` for differential testing against the oracle.
//! No silent failures: any error prints to stderr and exits non-zero.

use std::process::exit;

use clap::Parser;
use dockq_core::{score_pair, DockQError, InterfaceResult, RunOptions, RunResult};

#[derive(Parser, Debug)]
#[command(
    name = "dockq-rs",
    about = "DockQ — quality measure for protein/nucleic-acid docking models (Rust core)"
)]
struct Args {
    /// Path to the model file (PDB or mmCIF, optionally .gz).
    model: String,
    /// Path to the native/reference file (PDB or mmCIF, optionally .gz).
    native: String,

    /// Use the capri_peptide thresholds (CB interface, 4Å/8Å).
    #[arg(long = "capri_peptide")]
    capri_peptide: bool,
    /// Score a small-molecule pose (NOT supported in this build — hard-errors).
    #[arg(long = "small_molecule")]
    small_molecule: bool,
    /// Short single-line-per-interface output.
    #[arg(long)]
    short: bool,
    /// Verbose header.
    #[arg(long, short)]
    verbose: bool,
    /// Align by residue numbering instead of sequence.
    #[arg(long = "no_align")]
    no_align: bool,
    /// Number of worker threads (default: all cores).
    #[arg(long = "n_cpu")]
    n_cpu: Option<usize>,
    /// Accepted for compatibility; chunking is handled automatically by Rayon.
    #[arg(long = "max_chunk")]
    max_chunk: Option<usize>,
    /// Allowed sequence mismatches when clustering homologous chains.
    #[arg(long = "allowed_mismatches", default_value_t = 0)]
    allowed_mismatches: usize,
    /// Chain mapping MODELCHAINS:NATIVECHAINS (with optional `*` wildcards).
    #[arg(long)]
    mapping: Option<String>,
    /// Optimize the mapping on DockQ_F1 instead of DockQ (not implemented — hard-errors).
    #[arg(long = "optDockQF1")]
    opt_dockq_f1: bool,
    /// Write the full result as JSON to this file (like the reference `--json`).
    #[arg(long)]
    json: Option<String>,
    /// Print the oracle-compatible result JSON to stdout (for differential testing).
    #[arg(long = "diff_json")]
    diff_json: bool,
}

fn main() {
    let args = Args::parse();

    if args.opt_dockq_f1 {
        eprintln!(
            "dockq-rs error: --optDockQF1 (optimize mapping on DockQ_F1) is not implemented \
             in this build; it does not silently fall back to DockQ selection."
        );
        exit(2);
    }

    if let Some(n) = args.n_cpu {
        if n > 0 {
            // Best-effort; ignore error if a global pool already exists.
            let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
        }
    }

    let opts = RunOptions {
        no_align: args.no_align,
        capri_peptide: args.capri_peptide,
        small_molecule: args.small_molecule,
        mapping: args.mapping.clone(),
        allowed_mismatches: args.allowed_mismatches,
    };

    let result = match score_pair(&args.model, &args.native, &opts) {
        Ok(r) => r,
        Err(e) => {
            fail(&e);
        }
    };

    if let Some(path) = &args.json {
        let v = full_info_json(&result, args.capri_peptide, args.small_molecule);
        if let Err(e) = std::fs::write(path, serde_json::to_string(&v).unwrap()) {
            eprintln!("error writing JSON to {path}: {e}");
            exit(1);
        }
    }

    if args.diff_json {
        println!("{}", serde_json::to_string(&diff_json(&result)).unwrap());
        return;
    }

    print_results(&result, args.short, args.verbose, args.capri_peptide, args.small_molecule);
}

fn fail(e: &DockQError) -> ! {
    eprintln!("dockq-rs error: {e}");
    // Distinct exit code for the deferred small-molecule path.
    match e {
        DockQError::SmallMoleculeUnsupported => exit(3),
        _ => exit(1),
    }
}

fn score_name(capri_peptide: bool, small_molecule: bool) -> &'static str {
    if small_molecule {
        "DockQ-small_molecules"
    } else if capri_peptide {
        "DockQ-capri_peptide"
    } else {
        "DockQ"
    }
}

/// Reproduces the reference text output (header aside). For golden-file parity the header
/// lines all contain '*' and are stripped by the test harness's `grep -v "*"`.
fn print_results(r: &RunResult, short: bool, verbose: bool, capri_peptide: bool, small_molecule: bool) {
    let score = score_name(capri_peptide, small_molecule);
    let n = r.best_result.len();

    if short {
        println!(
            "Total {score} over {n} native interfaces: {:.3} with {} model:native mapping",
            r.global_dockq, r.best_mapping_str
        );
        for (chains, res) in &r.best_result {
            let c: Vec<char> = chains.chars().collect();
            let score_str = format!(
                "DockQ {:.3} iRMSD {:.3} LRMSD {:.3} fnat {:.3} fnonnat {:.3} F1 {:.3} clashes {}",
                res.dockq, res.irmsd, res.lrmsd, res.fnat, res.fnonnat, res.f1, res.clashes
            );
            println!(
                "{score_str} mapping {}{}:{}{} {} {} {} -> {} {} {}",
                res.chain1, res.chain2, c[0], c[1], r.model, res.chain1, res.chain2, r.native, c[0], c[1]
            );
        }
    } else {
        print_header(verbose, capri_peptide);
        println!("Model  : {}", r.model);
        println!("Native : {}", r.native);
        println!(
            "Total {score} over {n} native interfaces: {:.3} with {} model:native mapping",
            r.global_dockq, r.best_mapping_str
        );
        for (chains, res) in &r.best_result {
            let c: Vec<char> = chains.chars().collect();
            println!("Native chains: {}, {}", c[0], c[1]);
            println!("\tModel chains: {}, {}", res.chain1, res.chain2);
            println!("\tDockQ: {:.3}", res.dockq);
            println!("\tiRMSD: {:.3}", res.irmsd);
            println!("\tLRMSD: {:.3}", res.lrmsd);
            println!("\tfnat: {:.3}", res.fnat);
            println!("\tfnonnat: {:.3}", res.fnonnat);
            println!("\tF1: {:.3}", res.f1);
            println!("\tclashes: {}", res.clashes);
        }
    }
}

fn print_header(verbose: bool, capri_peptide: bool) {
    println!("****************************************************************");
    println!("*                       DockQ                                  *");
    println!("*   Scoring function for biomolecular docking models (Rust)     *");
    println!("*   DockQ score legend:                                         *");
    println!("*    0.00 <= DockQ <  0.23 - Incorrect                          *");
    println!("*    0.23 <= DockQ <  0.49 - Acceptable quality                 *");
    println!("*    0.49 <= DockQ <  0.80 - Medium quality                     *");
    println!("*            DockQ >= 0.80 - High quality                       *");
    if verbose {
        let contact = if capri_peptide { "4A" } else { "5A" };
        let iface = if capri_peptide { "8A (CB)" } else { "10A (all heavy atoms)" };
        println!("*   Contact <{contact} (Fnat), interface <{iface}                 *");
    }
    println!("****************************************************************");
}

fn iface_to_json(res: &InterfaceResult, native_pair: &str) -> serde_json::Value {
    let c: Vec<char> = native_pair.chars().collect();
    serde_json::json!({
        "DockQ": res.dockq,
        "F1": res.f1,
        "iRMSD": res.irmsd,
        "LRMSD": res.lrmsd,
        "fnat": res.fnat,
        "nat_correct": res.nat_correct,
        "nat_total": res.nat_total,
        "fnonnat": res.fnonnat,
        "nonnat_count": res.nonnat_count,
        "model_total": res.model_total,
        "clashes": res.clashes,
        "len1": res.len1,
        "len2": res.len2,
        "class1": res.class1,
        "class2": res.class2,
        "chain1": res.chain1,
        "chain2": res.chain2,
        "is_het": false,
        "native_chain1": c[0].to_string(),
        "native_chain2": c[1].to_string(),
    })
}

/// Oracle-compatible JSON (matches oracle_run.py): for differential testing.
fn diff_json(r: &RunResult) -> serde_json::Value {
    let mut best = serde_json::Map::new();
    for (pair, res) in &r.best_result {
        best.insert(pair.clone(), iface_to_json(res, pair));
    }
    serde_json::json!({
        "model": r.model,
        "native": r.native,
        "best_mapping_str": r.best_mapping_str,
        "best_dockq": r.best_dockq,
        "GlobalDockQ": r.global_dockq,
        "best_result": best,
    })
}

/// Full info JSON for the `--json FILE` option (mirrors the reference `info` dict).
fn full_info_json(r: &RunResult, _capri: bool, _smallmol: bool) -> serde_json::Value {
    let mut best = serde_json::Map::new();
    for (pair, res) in &r.best_result {
        best.insert(pair.clone(), iface_to_json(res, pair));
    }
    let mut mapping = serde_json::Map::new();
    for (native, model) in &r.best_mapping {
        mapping.insert(native.clone(), serde_json::Value::String(model.clone()));
    }
    serde_json::json!({
        "model": r.model,
        "native": r.native,
        "best_dockq": r.best_dockq,
        "GlobalDockQ": r.global_dockq,
        "best_mapping_str": r.best_mapping_str,
        "best_mapping": mapping,
        "best_result": best,
    })
}
