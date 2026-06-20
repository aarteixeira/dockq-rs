#!/usr/bin/env python
"""Ground-truth oracle: dump DockQ load_PDB parse results to JSON.

For every example structure file, load it via DockQ.DockQ.load_PDB and serialize
the full parsed chain/residue/atom tree, including float32-exact coordinates
(stored both as decimal and via float.hex) so a Rust port can be validated
byte-for-byte.

Runs each file twice:
  - small_molecule=False  -> oracle/dumps/parse/<filename>.json
  - small_molecule=True   -> oracle/dumps/parse_het/<filename>.json

Errors are caught per-file and recorded in the JSON "error" field (never abort).

Regenerate:
  source /Users/Andre.Teixeira/projects/DockQ/.venv-baseline/bin/activate
  python /Users/Andre.Teixeira/projects/dockq-rs/oracle/dump_parse.py
"""
import json
import os
import sys
import traceback

import numpy as np

from DockQ.DockQ import load_PDB

EXAMPLES_DIR = "/Users/Andre.Teixeira/projects/DockQ/examples"
OUT_BASE = "/Users/Andre.Teixeira/projects/dockq-rs/oracle/dumps"

# All example structure files named in the task.
EXAMPLE_FILES = [
    "1A2K_r_l_b.model.pdb",
    "1A2K_r_l_b.pdb",
    "1EXB_r_l_b.model.pdb",
    "1EXB_r_l_b.pdb",
    "1EXB.cif.gz",
    "dimer_dimer.model.pdb",
    "dimer_dimer.pdb",
    "model.pdb",
    "native.pdb",
    "6qwn-assembly1.cif.gz",
    "6qwn-assembly2.cif.gz",
    "1HHO_hem.cif",
    "2HHB_hem.cif",
]


def coord_to_hex(coord):
    """float.hex of each float32 component, preserving exact float32 bits.

    coord is a numpy float32 array [x, y, z]. float(c) widens the float32 to a
    Python float without changing its value, so float.hex captures the exact
    stored float32 magnitude.
    """
    return [float.hex(float(c)) for c in coord]


def dump_residue(res):
    """Serialize a single Biopython residue.

    res.id == (hetflag, resseq, icode). n_atoms counts atoms as returned by
    get_atoms(); n_unique_atom_ids mirrors DockQ.list_atoms_per_residue which
    uses set(a.id for a in res.get_unpacked_list()) (altloc-expanded).
    """
    het_flag, resseq, icode = res.id
    atoms = list(res.get_atoms())
    unpacked = res.get_unpacked_list()
    atom_records = []
    for atom in atoms:
        coord = atom.coord  # numpy float32 [x, y, z]
        atom_records.append(
            {
                "name": atom.id,
                "element": atom.element,
                "altloc": atom.get_altloc(),
                "coord": [float(c) for c in coord],
                "coord_hex": coord_to_hex(coord),
            }
        )
    return {
        "resname": res.resname,
        "het_flag": het_flag,
        "resseq": int(resseq),
        "icode": icode,
        "n_atoms": len(atoms),
        "n_unique_atom_ids": len(set(a.id for a in unpacked)),
        "atoms": atom_records,
    }


def dump_chain(chain):
    residues = [dump_residue(res) for res in chain]
    return {
        "id": chain.id,
        "sequence": chain.sequence,
        "is_het": chain.is_het,
        "n_residues": len(residues),
        "residues": residues,
    }


def dump_structure(path, small_molecule):
    """Load one file and return the full serializable dict.

    Errors are captured into the 'error' field; on error 'chains' is [].
    """
    record = {"path": path, "error": None, "chains": []}
    try:
        model = load_PDB(path, small_molecule=small_molecule)
        record["chains"] = [dump_chain(chain) for chain in model]
    except Exception:
        record["error"] = traceback.format_exc()
    return record


def main():
    os.makedirs(os.path.join(OUT_BASE, "parse"), exist_ok=True)
    os.makedirs(os.path.join(OUT_BASE, "parse_het"), exist_ok=True)

    errors = []
    for fname in EXAMPLE_FILES:
        path = os.path.join(EXAMPLES_DIR, fname)
        for small_molecule, subdir in ((False, "parse"), (True, "parse_het")):
            record = dump_structure(path, small_molecule)
            out_path = os.path.join(OUT_BASE, subdir, fname + ".json")
            with open(out_path, "w") as fp:
                json.dump(record, fp, sort_keys=True, indent=1)
            tag = f"{subdir}/{fname}"
            if record["error"] is not None:
                first_line = record["error"].strip().splitlines()[-1]
                errors.append((tag, first_line))
                print(f"[ERROR] {tag}: {first_line}", file=sys.stderr)
            else:
                nch = len(record["chains"])
                print(f"[ok]    {tag}: {nch} chains")

    print(f"\nWrote parse dumps for {len(EXAMPLE_FILES)} files x2 modes.")
    if errors:
        print(f"{len(errors)} (file,mode) combos errored:")
        for tag, msg in errors:
            print(f"  {tag}: {msg}")
    else:
        print("No errors.")


if __name__ == "__main__":
    main()
