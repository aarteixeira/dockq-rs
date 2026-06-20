//! dockq-rs CLI. Default mode is drop-in-compatible single-pair scoring (output matches
//! the reference DockQ, header aside). The `batch` subcommand scores many model/native
//! jobs in parallel (Rust). `--diff_json` emits oracle-compatible JSON for differential
//! testing. No silent failures: any error prints to stderr and exits non-zero; batch jobs
//! report per-job errors explicitly and the run exits non-zero if any failed.

use std::process::exit;

use clap::{Args, Parser, Subcommand, ValueEnum};
use dockq_core::{
    score_one_vs_many, score_pair, score_pairs, BatchOutcome, DockQError, InterfaceResult,
    RunOptions, RunResult,
};

#[derive(Parser, Debug)]
#[command(
    name = "dockq-rs",
    about = "DockQ — quality measure for protein/nucleic-acid docking models (Rust core)",
    // Default (no subcommand) parses the single-pair args below; `batch` is opt-in.
    args_conflicts_with_subcommands = true
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    #[command(flatten)]
    single: SingleArgs,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Score many model/native jobs in parallel (one native vs many models, or pair lists).
    Batch(BatchArgs),
}

#[derive(Args, Debug)]
struct SingleArgs {
    /// Path to the model file (PDB or mmCIF, optionally .gz). Required unless using `batch`.
    model: Option<String>,
    /// Path to the native/reference file (PDB or mmCIF, optionally .gz).
    native: Option<String>,

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

#[derive(Args, Debug)]
struct BatchArgs {
    /// Native/reference structure to score every model against (one-vs-many modes).
    #[arg(long)]
    native: Option<String>,
    /// Model files to score against --native (space-separated, repeatable).
    #[arg(long, num_args = 1..)]
    models: Vec<String>,
    /// Directory: score every .pdb/.cif (optionally .gz) file in it against --native.
    #[arg(long = "models_dir")]
    models_dir: Option<String>,
    /// File listing model paths, one per line, to score against --native.
    #[arg(long = "models_from")]
    models_from: Option<String>,
    /// File listing "model native" pairs (whitespace-separated), one per line.
    #[arg(long = "pairs_from")]
    pairs_from: Option<String>,

    /// Output file (default: stdout).
    #[arg(long, short)]
    output: Option<String>,
    /// Output format.
    #[arg(long, value_enum, default_value_t = BatchFormat::Tsv)]
    format: BatchFormat,
    /// Sort rows by GlobalDockQ descending (handy for model ranking).
    #[arg(long)]
    sort: bool,

    // Scoring flags applied to every job.
    #[arg(long = "capri_peptide")]
    capri_peptide: bool,
    #[arg(long = "no_align")]
    no_align: bool,
    #[arg(long = "allowed_mismatches", default_value_t = 0)]
    allowed_mismatches: usize,
    /// Chain mapping applied to every job (with optional `*` wildcards).
    #[arg(long)]
    mapping: Option<String>,
    /// Number of worker threads (default: all cores).
    #[arg(long = "n_cpu")]
    n_cpu: Option<usize>,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum BatchFormat {
    Tsv,
    Json,
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Some(Commands::Batch(args)) => run_batch(args),
        None => run_single(cli.single),
    }
}

// ---------------------------------------------------------------------------
// Single-pair mode
// ---------------------------------------------------------------------------

fn run_single(args: SingleArgs) {
    if args.opt_dockq_f1 {
        eprintln!(
            "dockq-rs error: --optDockQF1 (optimize mapping on DockQ_F1) is not implemented \
             in this build; it does not silently fall back to DockQ selection."
        );
        exit(2);
    }
    set_threads(args.n_cpu);
    let _ = args.max_chunk;

    let (model, native) = match (&args.model, &args.native) {
        (Some(m), Some(n)) => (m, n),
        _ => {
            eprintln!(
                "dockq-rs error: a <model> and <native> are required.\n\
                 Usage: dockq-rs <model> <native> [OPTIONS]   (or: dockq-rs batch ...)"
            );
            exit(2);
        }
    };

    let opts = RunOptions {
        no_align: args.no_align,
        capri_peptide: args.capri_peptide,
        small_molecule: args.small_molecule,
        mapping: args.mapping.clone(),
        allowed_mismatches: args.allowed_mismatches,
    };

    let result = match score_pair(model, native, &opts) {
        Ok(r) => r,
        Err(e) => fail(&e),
    };

    if let Some(path) = &args.json {
        let v = full_info_json(&result);
        if let Err(e) = std::fs::write(path, serde_json::to_string(&v).unwrap()) {
            eprintln!("error writing JSON to {path}: {e}");
            exit(1);
        }
    }

    if args.diff_json {
        println!("{}", serde_json::to_string(&diff_json(&result)).unwrap());
        return;
    }

    print_results(
        &result,
        args.short,
        args.verbose,
        args.capri_peptide,
        args.small_molecule,
    );
}

// ---------------------------------------------------------------------------
// Batch mode
// ---------------------------------------------------------------------------

fn run_batch(args: BatchArgs) {
    set_threads(args.n_cpu);
    let opts = RunOptions {
        no_align: args.no_align,
        capri_peptide: args.capri_peptide,
        small_molecule: false,
        mapping: args.mapping.clone(),
        allowed_mismatches: args.allowed_mismatches,
    };

    let outcomes: Vec<BatchOutcome> = if let Some(pairs_file) = &args.pairs_from {
        if args.native.is_some() || !args.models.is_empty() || args.models_dir.is_some() {
            eprintln!("dockq-rs batch error: --pairs_from is mutually exclusive with --native/--models/--models_dir");
            exit(2);
        }
        let pairs = read_pairs(pairs_file);
        score_pairs(&pairs, &opts)
    } else {
        let native = match &args.native {
            Some(n) => n,
            None => {
                eprintln!(
                    "dockq-rs batch error: provide --native with --models/--models_dir/--models_from, \
                     or use --pairs_from"
                );
                exit(2);
            }
        };
        let mut models: Vec<String> = args.models.clone();
        if let Some(dir) = &args.models_dir {
            match dockq_core::scan_structures(dir) {
                Ok(found) => models.extend(found),
                Err(e) => {
                    eprintln!("dockq-rs batch error: scanning {dir}: {e}");
                    exit(2);
                }
            }
        }
        if let Some(f) = &args.models_from {
            models.extend(read_lines(f));
        }
        if models.is_empty() {
            eprintln!("dockq-rs batch error: no models (use --models / --models_dir / --models_from)");
            exit(2);
        }
        score_one_vs_many(native, &models, &opts)
    };

    // Display order (input order, or sorted by GlobalDockQ descending).
    let mut order: Vec<usize> = (0..outcomes.len()).collect();
    if args.sort {
        order.sort_by(|&a, &b| {
            global_of(&outcomes[b])
                .partial_cmp(&global_of(&outcomes[a]))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    let text = match args.format {
        BatchFormat::Tsv => batch_tsv(&outcomes, &order),
        BatchFormat::Json => batch_json(&outcomes, &order),
    };
    match &args.output {
        Some(path) => {
            if let Err(e) = std::fs::write(path, text) {
                eprintln!("error writing output to {path}: {e}");
                exit(1);
            }
        }
        None => print!("{text}"),
    }

    let n_err = outcomes.iter().filter(|o| o.result.is_err()).count();
    if n_err > 0 {
        eprintln!("dockq-rs batch: {n_err}/{} job(s) failed", outcomes.len());
        exit(1);
    }
}

/// GlobalDockQ for sorting (failed jobs sink to the bottom).
fn global_of(o: &BatchOutcome) -> f64 {
    match &o.result {
        Ok(r) => r.global_dockq,
        Err(_) => f64::NEG_INFINITY,
    }
}

fn batch_tsv(outcomes: &[BatchOutcome], order: &[usize]) -> String {
    let mut s = String::from("model\tnative\tGlobalDockQ\tDockQ_sum\tn_interfaces\tmapping\tstatus\n");
    for &i in order {
        let o = &outcomes[i];
        match &o.result {
            Ok(r) => s.push_str(&format!(
                "{}\t{}\t{:.6}\t{:.6}\t{}\t{}\tok\n",
                o.model,
                o.native,
                r.global_dockq,
                r.best_dockq,
                r.best_result.len(),
                r.best_mapping_str
            )),
            Err(e) => s.push_str(&format!(
                "{}\t{}\t\t\t\t\terror: {}\n",
                o.model,
                o.native,
                e.to_string().replace(['\t', '\n'], " ")
            )),
        }
    }
    s
}

fn batch_json(outcomes: &[BatchOutcome], order: &[usize]) -> String {
    let arr: Vec<serde_json::Value> = order
        .iter()
        .map(|&i| {
            let o = &outcomes[i];
            match &o.result {
                Ok(r) => serde_json::json!({
                    "model": o.model, "native": o.native, "ok": true,
                    "GlobalDockQ": r.global_dockq, "best_dockq": r.best_dockq,
                    "best_mapping_str": r.best_mapping_str,
                    "best_result": result_map_json(r),
                }),
                Err(e) => serde_json::json!({
                    "model": o.model, "native": o.native, "ok": false,
                    "error": e.to_string(),
                }),
            }
        })
        .collect();
    serde_json::to_string_pretty(&serde_json::Value::Array(arr)).unwrap()
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn set_threads(n_cpu: Option<usize>) {
    if let Some(n) = n_cpu {
        if n > 0 {
            // Best-effort; ignore error if a global pool already exists.
            let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
        }
    }
}

/// Read "model native" pairs (whitespace-separated), skipping blanks and `#` comments.
fn read_pairs(path: &str) -> Vec<(String, String)> {
    let text = read_file_or_die(path);
    let mut pairs = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let toks: Vec<&str> = t.split_whitespace().collect();
        if toks.len() != 2 {
            eprintln!(
                "dockq-rs batch error: {path}:{}: expected 'model native', got {:?}",
                lineno + 1,
                t
            );
            exit(2);
        }
        pairs.push((toks[0].to_string(), toks[1].to_string()));
    }
    pairs
}

/// Read non-blank, non-comment lines as paths.
fn read_lines(path: &str) -> Vec<String> {
    read_file_or_die(path)
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(String::from)
        .collect()
}

fn read_file_or_die(path: &str) -> String {
    match std::fs::read_to_string(path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("dockq-rs batch error: reading {path}: {e}");
            exit(2);
        }
    }
}

fn fail(e: &DockQError) -> ! {
    eprintln!("dockq-rs error: {e}");
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
fn print_results(
    r: &RunResult,
    short: bool,
    verbose: bool,
    capri_peptide: bool,
    small_molecule: bool,
) {
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

fn result_map_json(r: &RunResult) -> serde_json::Map<String, serde_json::Value> {
    let mut best = serde_json::Map::new();
    for (pair, res) in &r.best_result {
        best.insert(pair.clone(), iface_to_json(res, pair));
    }
    best
}

/// Oracle-compatible JSON (matches oracle_run.py): for differential testing.
fn diff_json(r: &RunResult) -> serde_json::Value {
    serde_json::json!({
        "model": r.model,
        "native": r.native,
        "best_mapping_str": r.best_mapping_str,
        "best_dockq": r.best_dockq,
        "GlobalDockQ": r.global_dockq,
        "best_result": result_map_json(r),
    })
}

/// Full info JSON for the `--json FILE` option (mirrors the reference `info` dict).
fn full_info_json(r: &RunResult) -> serde_json::Value {
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
        "best_result": result_map_json(r),
    })
}
