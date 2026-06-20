#!/usr/bin/env python
"""Ground-truth oracle: dump DockQ per-interface intermediates to JSON.

For a curated list of (model, native) cases, this replicates DockQ.main()'s
chain-mapping search to obtain the BEST chain_map, then for that map dumps, for
each scored native chain pair:

  - alignments      : {seqA, matches, seqB} from format_alignment(align_chains(...))
  - residue_distances: shape, contact counts at the relevant squared thresholds,
                       sha256 of the float32 matrix .tobytes(), and the full
                       nested list when small enough
  - fnat_stats      : [n_shared, n_nonnat, n_native, n_model] from get_fnat_stats
  - irms, lrms      : floats (recomputed via the same path calc_DockQ uses)
  - info            : the full info dict returned by run_on_all_native_interfaces

All numeric library functions are reused (no reimplementation of DockQ logic).

Regenerate (with the reference DockQ on PYTHONPATH; set DOCKQ_REPO if it is not a
sibling ../DockQ):
  python oracle/dump_interfaces.py
"""
import hashlib
import itertools
import json
import os
import sys
import traceback

import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer

from DockQ.DockQ import (
    align_chains,
    calc_DockQ,
    count_chain_combinations,
    format_alignment,
    format_mapping,
    format_mapping_string,
    get_aligned_residues,
    get_all_chain_maps,
    get_interacting_pairs,
    get_residue_distances,
    group_chains,
    load_PDB,
    run_on_all_native_interfaces,
    subset_atoms,
)
from DockQ.operations import get_fnat_stats
from DockQ.constants import (
    BACKBONE_ATOMS,
    FNAT_THRESHOLD,
    FNAT_THRESHOLD_PEPTIDE,
    INTERFACE_THRESHOLD,
    INTERFACE_THRESHOLD_PEPTIDE,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_DOCKQ = os.environ.get("DOCKQ_REPO") or os.path.join(
    os.path.dirname(os.path.dirname(_HERE)), "DockQ"
)
EX = os.environ.get("DOCKQ_EXAMPLES", os.path.join(_DOCKQ, "examples"))
OUT = os.path.join(_HERE, "dumps", "interfaces")

# Full-matrix dump cap: emit nested list only when r*c <= this.
FULL_MATRIX_CAP = 4000

# case name -> (model_path, native_path, kwargs)
# kwargs may contain: no_align, capri_peptide, allowed_mismatches, mapping
CASES = {
    "1A2K": (f"{EX}/1A2K_r_l_b.model.pdb", f"{EX}/1A2K_r_l_b.pdb", {}),
    "1A2K_noalign": (
        f"{EX}/1A2K_r_l_b.model.pdb",
        f"{EX}/1A2K_r_l_b.pdb",
        {"no_align": True},
    ),
    "dimer_dimer": (f"{EX}/dimer_dimer.model.pdb", f"{EX}/dimer_dimer.pdb", {}),
    # model.pdb/native.pdb have slightly different sequences; CLI uses
    # --allowed_mismatches 1 (see run_test.sh line 25).
    "model_native": (
        f"{EX}/model.pdb",
        f"{EX}/native.pdb",
        {"allowed_mismatches": 1},
    ),
    "1EXB": (f"{EX}/1EXB_r_l_b.model.pdb", f"{EX}/1EXB_r_l_b.pdb", {}),
    "1EXB_cif": (f"{EX}/1EXB_r_l_b.model.pdb", f"{EX}/1EXB.cif.gz", {}),
    # CLI: 6qwn-assembly1 (model) vs 6qwn-assembly2 (native) --capri_peptide
    # (run_test.sh line 47; golden testdata/6q2n_peptide.dockq).
    "6qwn_capri": (
        f"{EX}/6qwn-assembly1.cif.gz",
        f"{EX}/6qwn-assembly2.cif.gz",
        {"capri_peptide": True},
    ),
}


def matrix_to_hex_sha(mat_f32):
    """sha256 of the row-major float32 matrix bytes (mat must already be f32)."""
    return hashlib.sha256(mat_f32.tobytes()).hexdigest()


def dump_residue_distances(dist, fnat_threshold, interface_threshold):
    """Serialize a residue-distance matrix (squared distances, Angstrom^2).

    Counts use squared thresholds, matching DockQ's <threshold**2 comparisons.
    n_lt_25/100/4 are computed at the fixed 5/10/2 A reference (squared 25/100/4)
    for stable cross-checking; n_lt_fnat/n_lt_interface use the case's actual
    thresholds (which differ for capri_peptide).
    """
    d = np.asarray(dist)
    d_f32 = d.astype(np.float32)
    r, c = d.shape
    full = d.tolist() if (r * c <= FULL_MATRIX_CAP) else None
    return {
        "shape": [int(r), int(c)],
        "n_lt_25": int(np.count_nonzero(d < 25.0)),  # < 5^2
        "n_lt_100": int(np.count_nonzero(d < 100.0)),  # < 10^2
        "n_lt_4": int(np.count_nonzero(d < 4.0)),  # < 2^2 (clash)
        "n_lt_fnat": int(np.count_nonzero(d < fnat_threshold ** 2)),
        "n_lt_interface": int(np.count_nonzero(d < interface_threshold ** 2)),
        "fnat_threshold": float(fnat_threshold),
        "interface_threshold": float(interface_threshold),
        "dtype": str(d.dtype),
        "bytes_sha256": matrix_to_hex_sha(d_f32),
        "full": full,
    }


def recompute_protein_interface(model_chains, native_chains, capri_peptide):
    """Recompute, via DockQ's own helpers, the same intermediates calc_DockQ
    uses for a protein-protein native pair: alignments, the sample/native
    residue-distance matrix that feeds get_fnat_stats, fnat_stats, irms, lrms.

    This mirrors calc_DockQ (DockQ.py) step-for-step using the library funcs so
    the dumped values are exactly what DockQ scores on.
    Returns a dict, or None if the native pair has no interface (nat_total==0).
    """
    fnat_threshold = FNAT_THRESHOLD if not capri_peptide else FNAT_THRESHOLD_PEPTIDE
    interface_threshold = (
        INTERFACE_THRESHOLD if not capri_peptide else INTERFACE_THRESHOLD_PEPTIDE
    )

    # Alignments (same as run_on_chains): use_numbering is False here because
    # all curated protein cases run with no_align=False except 1A2K_noalign,
    # which is handled by passing use_numbering through the caller.
    alignments = []
    for mc, nc in zip(model_chains, native_chains):
        aln = align_chains(mc, nc, use_numbering=recompute_protein_interface.no_align)
        alignments.append(tuple(format_alignment(aln).values()))
    alignments = tuple(alignments)

    # native contacts on untouched native (calc_DockQ).
    ref_res_distances = get_residue_distances(native_chains[0], native_chains[1], "ref")
    nat_total = np.nonzero(np.asarray(ref_res_distances) < fnat_threshold ** 2)[0].shape[0]
    if nat_total == 0:
        return None

    aligned_sample_1, aligned_ref_1 = get_aligned_residues(
        model_chains[0], native_chains[0], alignments[0]
    )
    aligned_sample_2, aligned_ref_2 = get_aligned_residues(
        model_chains[1], native_chains[1], alignments[1]
    )

    sample_res_distances = get_residue_distances(
        aligned_sample_1, aligned_sample_2, "sample"
    )
    if ref_res_distances.shape != sample_res_distances.shape:
        ref_res_distances = get_residue_distances(aligned_ref_1, aligned_ref_2, "ref")

    nat_correct, nonnat_count, n_native, model_total = get_fnat_stats(
        sample_res_distances, ref_res_distances, threshold=fnat_threshold
    )

    # iRMSD: interface defined on (capri-adjusted) native distances.
    ref_for_interface = ref_res_distances
    if capri_peptide:
        ref_for_interface = get_residue_distances(
            aligned_ref_1, aligned_ref_2, "ref", all_atom=False
        )
    interacting_pairs = get_interacting_pairs(
        ref_for_interface, threshold=interface_threshold ** 2
    )
    s_int1, r_int1 = subset_atoms(
        aligned_sample_1, aligned_ref_1, atom_types=BACKBONE_ATOMS,
        residue_subset=interacting_pairs[0],
    )
    s_int2, r_int2 = subset_atoms(
        aligned_sample_2, aligned_ref_2, atom_types=BACKBONE_ATOMS,
        residue_subset=interacting_pairs[1],
    )
    sample_interface_atoms = np.asarray(s_int1 + s_int2)
    ref_interface_atoms = np.asarray(r_int1 + r_int2)
    si = SVDSuperimposer()
    si.set(sample_interface_atoms, ref_interface_atoms)
    si.run()
    irms = si.get_rms()

    # LRMSD: superpose on receptor backbone, rms on ligand backbone.
    ref_group1_size = len(native_chains[0])
    ref_group2_size = len(native_chains[1])
    receptor_chains = (
        (aligned_ref_1, aligned_sample_1)
        if ref_group1_size > ref_group2_size
        else (aligned_ref_2, aligned_sample_2)
    )
    ligand_chains = (
        (aligned_ref_1, aligned_sample_1)
        if ref_group1_size <= ref_group2_size
        else (aligned_ref_2, aligned_sample_2)
    )
    rec_native, rec_sample = subset_atoms(
        receptor_chains[0], receptor_chains[1], atom_types=BACKBONE_ATOMS, what="receptor"
    )
    lig_native, lig_sample = subset_atoms(
        ligand_chains[0], ligand_chains[1], atom_types=BACKBONE_ATOMS, what="ligand"
    )
    si2 = SVDSuperimposer()
    si2.set(np.asarray(rec_native), np.asarray(rec_sample))
    si2.run()
    rot, tran = si2.get_rotran()
    rotated = np.dot(np.asarray(lig_sample), rot) + tran
    lrms = si2._rms(np.asarray(lig_native), rotated)

    return {
        "alignments": [
            {"seqA": a[0], "matches": a[1], "seqB": a[2]} for a in alignments
        ],
        "residue_distances": dump_residue_distances(
            sample_res_distances, fnat_threshold, interface_threshold
        ),
        "native_residue_distances": dump_residue_distances(
            ref_res_distances, fnat_threshold, interface_threshold
        ),
        "fnat_stats": [
            int(nat_correct), int(nonnat_count), int(n_native), int(model_total)
        ],
        "irms": float(irms),
        "lrms": float(lrms),
    }


recompute_protein_interface.no_align = False


def best_chain_map(model_structure, native_structure, kwargs):
    """Replicate DockQ.main()'s chain-mapping search; return (best_map, best_dockq).

    Serial deterministic loop (n_cpu=1 equivalent): iterates the same chain_maps
    iterator main() builds and selects by total_dockq > best (first-wins ties),
    identical to main()'s reduction over the parallel results.
    """
    mapping_str = kwargs.get("mapping")
    small_molecule = kwargs.get("small_molecule", False)
    allowed_mismatches = kwargs.get("allowed_mismatches", 0)
    no_align = kwargs.get("no_align", False)
    capri_peptide = kwargs.get("capri_peptide", False)

    initial_mapping, model_chains, native_chains = format_mapping(
        mapping_str, small_molecule
    )
    model_chains = (
        [c.id for c in model_structure] if not model_chains else model_chains
    )
    native_chains = (
        [c.id for c in native_structure] if not native_chains else native_chains
    )

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
        allowed_mismatches,
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
    best_mapping = None
    if num_chain_combinations > 1:
        for chain_map in chain_maps_:
            _, total_dockq = run_on_all_native_interfaces(
                model_structure,
                native_structure,
                chain_map=chain_map,
                no_align=no_align,
                capri_peptide=capri_peptide,
                low_memory=low_memory,
            )
            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_mapping = chain_map
    else:
        best_mapping = next(chain_maps)
        _, best_dockq = run_on_all_native_interfaces(
            model_structure,
            native_structure,
            chain_map=best_mapping,
            no_align=no_align,
            capri_peptide=capri_peptide,
            low_memory=False,
        )
    return best_mapping, best_dockq


def dump_case(name, model_path, native_path, kwargs):
    record = {
        "name": name,
        "model": model_path,
        "native": native_path,
        "kwargs": kwargs,
        "error": None,
    }
    try:
        small_molecule = kwargs.get("small_molecule", False)
        no_align = kwargs.get("no_align", False)
        capri_peptide = kwargs.get("capri_peptide", False)

        model_structure = load_PDB(
            model_path, chains=[], small_molecule=small_molecule
        )
        native_structure = load_PDB(
            native_path, chains=[], small_molecule=small_molecule
        )

        best_map, best_dockq = best_chain_map(
            model_structure, native_structure, kwargs
        )
        record["best_chain_map"] = best_map
        record["best_mapping_str"] = format_mapping_string(best_map)
        record["best_total_dockq"] = float(best_dockq)

        # Full info via the library (authoritative).
        result_mapping, total_dockq = run_on_all_native_interfaces(
            model_structure,
            native_structure,
            chain_map=best_map,
            no_align=no_align,
            capri_peptide=capri_peptide,
            low_memory=False,
        )
        record["total_dockq"] = float(total_dockq)
        record["n_interfaces"] = len(result_mapping)
        record["GlobalDockQ"] = (
            float(total_dockq) / len(result_mapping) if result_mapping else None
        )

        # Per native-pair dump.
        recompute_protein_interface.no_align = no_align
        interfaces = {}
        native_chain_ids = list(best_map.keys())
        for chain_pair in itertools.combinations(native_chain_ids, 2):
            key = "".join(chain_pair)
            if key not in result_mapping:
                continue  # null interface or skipped (e.g. same model chain)
            info = result_mapping[key]
            entry = {"info": json_safe(info)}

            native_chains = tuple(native_structure[c] for c in chain_pair)
            model_chains = tuple(
                model_structure[best_map[c]] for c in chain_pair
            )

            is_het = bool(native_chains[0].is_het or native_chains[1].is_het)
            entry["is_het"] = is_het
            if not is_het:
                detail = recompute_protein_interface(
                    model_chains, native_chains, capri_peptide
                )
                if detail is not None:
                    entry.update(detail)
            else:
                # Small-molecule path: info already carries DockQ/LRMSD/mapping.
                # Dump the alignments used (receptor alignment is meaningful).
                alignments = []
                for mc, nc in zip(model_chains, native_chains):
                    aln = align_chains(mc, nc, use_numbering=no_align)
                    alignments.append(format_alignment(aln))
                entry["alignments"] = [
                    {"seqA": a["seqA"], "matches": a["matches"], "seqB": a["seqB"]}
                    for a in alignments
                ]
                entry["lrms"] = float(info.get("LRMSD")) if info.get("LRMSD") is not None else None
            interfaces[key] = entry

        record["interfaces"] = interfaces
    except Exception:
        record["error"] = traceback.format_exc()
    return record


def json_safe(obj):
    """Make an info dict JSON-serializable: drop the chain_map diagnostic dict's
    object refs are fine (they are plain str:str), and cast numpy scalars.
    """
    out = {}
    for k, v in obj.items():
        if isinstance(v, (np.integer,)):
            out[k] = int(v)
        elif isinstance(v, (np.floating,)):
            out[k] = float(v)
        elif isinstance(v, np.ndarray):
            out[k] = v.tolist()
        elif isinstance(v, dict):
            out[k] = {str(kk): (str(vv) if not isinstance(vv, (int, float, str, type(None))) else vv) for kk, vv in v.items()}
        else:
            out[k] = v
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    errors = []
    for name, (model_path, native_path, kwargs) in CASES.items():
        record = dump_case(name, model_path, native_path, kwargs)
        out_path = os.path.join(OUT, f"{name}.json")
        with open(out_path, "w") as fp:
            json.dump(record, fp, sort_keys=True, indent=1)
        if record["error"] is not None:
            last = record["error"].strip().splitlines()[-1]
            errors.append((name, last))
            print(f"[ERROR] {name}: {last}", file=sys.stderr)
        else:
            print(
                f"[ok]    {name}: best={record['best_mapping_str']} "
                f"n_iface={record['n_interfaces']} "
                f"GlobalDockQ={record['GlobalDockQ']:.6f}"
            )

    print(f"\nWrote interface dumps for {len(CASES)} cases.")
    if errors:
        print(f"{len(errors)} cases errored:")
        for name, msg in errors:
            print(f"  {name}: {msg}")
    else:
        print("No errors.")


if __name__ == "__main__":
    main()
