"""Verify frame times and pocket distances from a production trajectory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array

from reproduce import read_rows


DATA_DIR = Path(__file__).resolve().parents[1] / "assets"
DISTANCE_COLUMN = "fixed_pocket_heavy_minimum_image_distance_nm"
EPIMER_STEM = "D-Ru5P_C4-epimer"
DISTANCE_ATOL_NM = 1e-6
GROMACS_DISTANCE_ATOL_NM = 5e-6


def trajectory_files(args, parser, systems):
    if args.system_dir is not None:
        if args.system is not None:
            parser.error("--system applies to --trajectory-dir")
        return args.system_dir.name, args.system_dir / "md.tpr", args.system_dir / "md.xtc"

    if args.system is None:
        parser.error("--trajectory-dir requires --system")
    if args.system not in systems:
        parser.error(f"Choose a system from {', '.join(systems)}")
    system = args.system
    zenodo_stem = EPIMER_STEM if system == "substrate" else system
    path = args.trajectory_dir / f"RPE_{zenodo_stem}_100ns"
    return system, path.with_name(path.name + ".tpr"), path.with_name(path.name + ".xtc")


def archived_atom_groups(raw, mapping_file, pocket):
    mapping = read_rows(mapping_file)
    indices = [int(row["source_atom"]) - 1 for row in mapping]
    mapped = raw.atoms[indices]
    if not np.array_equal(mapped.names, [row["name"] for row in mapping]):
        raise ValueError("Production topology atom names differ from the reference mapping")
    if not np.array_equal(mapped.resnames, [row["residue_name"] for row in mapping]):
        raise ValueError("Production topology residue names differ from the reference mapping")
    protein = np.array([row["chain"] == "A" for row in mapping])
    heavy = np.array([row["element"] != "H" for row in mapping])
    in_pocket = np.isin([int(row["residue_number"]) for row in mapping], pocket)
    return mapped[protein & heavy & in_pocket], mapped[~protein & heavy]


def native_atom_groups(raw, pocket):
    protein = raw.residues[:282].atoms
    ligand = raw.residues[282:283].atoms
    pocket_atoms = protein[np.isin(protein.resids, pocket) & ~np.char.startswith(protein.names.astype(str), "H")]
    ligand_atoms = ligand[~np.char.startswith(ligand.names.astype(str), "H")]
    if set(pocket_atoms.resids) != set(pocket) or len(ligand_atoms) != 14 or set(ligand.resnames) != {"MOL"}:
        raise ValueError("Production topology does not match the fixed pocket and native substrate")
    return pocket_atoms, ligand_atoms


def verify(data, system, topology, trajectory):
    track = data / "md_observables"
    pocket = [int(row["resid"]) for row in read_rows(data / "docking/fixed_pocket_definition.tsv")]
    expected = read_rows(track / f"md_metrics_{system}.tsv")
    raw = mda.Universe(str(topology), str(trajectory))
    native = system.startswith("native_D-Ru5P_")
    try:
        if native:
            pocket_atoms, ligand_atoms = native_atom_groups(raw, pocket)
        else:
            pocket_atoms, ligand_atoms = archived_atom_groups(raw, track / "atom_maps" / f"{system}.tsv", pocket)

        times = []
        distances = []
        for frame in raw.trajectory:
            times.append(frame.time / 1000)
            if len(ligand_atoms):
                values = distance_array(pocket_atoms.positions, ligand_atoms.positions, box=frame.dimensions)
                distances.append(values.min() / 10)
    finally:
        raw.trajectory.close()

    expected_times = [float(row["time_ns"]) for row in expected]
    np.testing.assert_allclose(times, expected_times, rtol=0, atol=1e-6)
    summary = {"system": system, "frames": len(times)}
    if distances:
        arr = np.asarray(distances)
        expected_distances = np.array([float(row[DISTANCE_COLUMN]) for row in expected])
        # The paired reference series was exported by GROMACS in nanometers.
        tolerance = GROMACS_DISTANCE_ATOL_NM if native else DISTANCE_ATOL_NM
        np.testing.assert_allclose(arr, expected_distances, rtol=0, atol=tolerance)
        window = arr[(np.asarray(times) >= 50) & (np.asarray(times) <= 100)]
        summary.update(
            max_distance_difference_nm=float(np.max(np.abs(arr - expected_distances))),
            frames_50_to_100_ns=len(window),
            frames_below_4_angstrom=int(np.sum(window < 0.4)),
            frames_below_5_angstrom=int(np.sum(window < 0.5)),
        )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA_DIR, help="Directory of analysis tables")
    location = parser.add_mutually_exclusive_group(required=True)
    location.add_argument("--trajectory-dir", type=Path, help="Directory of downloaded Zenodo files")
    location.add_argument("--system-dir", type=Path, help="Directory containing md.tpr and md.xtc")
    parser.add_argument("--system", help="System identifier used in the reference tables")
    args = parser.parse_args()
    systems = sorted(path.stem.removeprefix("md_metrics_") for path in (args.data / "md_observables").glob("md_metrics_*.tsv"))
    system, topology, trajectory = trajectory_files(args, parser, systems)
    print(json.dumps(verify(args.data, system, topology, trajectory), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
