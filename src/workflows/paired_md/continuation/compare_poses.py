"""Check the final Ru5P configurations and compare core-aligned ligand poses."""

import argparse
import csv
import math
from pathlib import Path
from statistics import fmean, median
import subprocess


def gro_atoms(path):
    with path.open(encoding="ascii") as source:
        next(source)
        count = int(next(source))
        return [(int(line[:5]), line[5:10].strip(), line[10:15].strip(),
                 tuple(float(line[i:i+8]) for i in (20, 28, 36)))
                for _, line in zip(range(count), source)]


def neighbor_volume(coords, name, neighbors):
    x, y, z = [tuple(a-b for a, b in zip(coords[atom], coords[name])) for atom in neighbors]
    return sum(x[i] * (y[(i+1) % 3] * z[(i+2) % 3] - y[(i+2) % 3] * z[(i+1) % 3])
               for i in range(3))


def chiral_volumes(atoms):
    ligand = {name: coord for _, residue, name, coord in atoms if residue == "MOL"}
    assert len(ligand) == 23
    return (neighbor_volume(ligand, "C3", ("O3", "C2", "C4")),
            neighbor_volume(ligand, "C4", ("O4", "C3", "C5")))


def indices(path, groups):
    with path.open("w", encoding="ascii") as out:
        for name, atoms in groups.items():
            out.write(f"[ {name} ]\n")
            for start in range(0, len(atoms), 15):
                out.write(" ".join(map(str, atoms[start:start+15])) + "\n")


def assess(gmx, run, output, spacing, duration):
    start = gro_atoms(run / "setup/start.gro")
    final = gro_atoms(run / "pilot/pilot.gro")
    before = chiral_volumes(start)
    after = chiral_volumes(final)
    if not all(a < 0 and b < 0 for a, b in zip(before, after)):
        raise ValueError(f"Unexpected Ru5P tetrahedral geometry in {run}")
    core = [i for i, (resid, resname, name, _) in enumerate(start, 1)
            if 48 <= resid <= 279 and name == "CA" and resname != "MOL"]
    ligand = [i for i, (_, resname, name, _) in enumerate(start, 1)
              if resname == "MOL" and not name.lstrip("0123456789").startswith("H")]
    if len(core) != 232 or len(ligand) != 14:
        raise ValueError(f"Unexpected core or Ru5P atom selection in {run}")
    output.mkdir()
    solute = [i for i, (_, resname, _, _) in enumerate(start, 1) if resname not in {"SOL", "NA", "CL"}]
    if solute != list(range(1, len(solute) + 1)):
        raise ValueError(f"Solute atom numbering is not continuous in {run}")
    imaging_groups = output / "imaging.ndx"
    indices(imaging_groups, {"CoreCa": core, "Solute": solute})
    imaged = output / "solute_centered.xtc"
    reference = output / "initial_centered.gro"
    with (output / "imaging.log").open("w", encoding="utf-8") as log:
        command = [gmx, "trjconv", "-s", str(run / "pilot/pilot.tpr"),
                   "-f", str(run / "pilot/pilot.xtc"), "-n", str(imaging_groups),
                   "-pbc", "mol", "-center", "-ur", "compact"]
        for target, options in ((imaged, []), (reference, ["-dump", "0"])):
            subprocess.run(command + ["-o", str(target)] + options, input="0\n1\n", text=True,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
    groups = output / "rmsd.ndx"
    series = output / "ligand_core_aligned_rmsd_nm.xvg"
    indices(groups, {"CoreCa": core, "Ru5PHeavy": ligand})
    with (output / "rmsd.log").open("w", encoding="utf-8") as log:
        subprocess.run([gmx, "rms", "-s", str(reference), "-f", str(imaged), "-n", str(groups),
                        "-o", str(series), "-tu", "ps", "-xvg", "none", "-nomw"],
                       input="0\n1\n", text=True, stdout=log, stderr=subprocess.STDOUT, check=True)
    rows = [tuple(map(float, line.split())) for line in series.read_text().splitlines()
            if line.strip() and not line.startswith(("#", "@"))]
    times = range(0, duration + 1, spacing)
    if len(rows) != len(times) or any(abs(t-ref) > 0.05 or not math.isfinite(rms) or rms < 0
                                     for (t, rms), ref in zip(rows, times)):
        raise ValueError(f"Unexpected RMSD series in {series}")
    rmsd = [distance for _, distance in rows]
    if rmsd[0] > 0.005:
        raise ValueError(f"RMSD reference frame does not match the first imaged frame in {series}")
    late = [distance for time, distance in rows if time >= duration / 2]
    return (len(rows), rmsd[-1], fmean(rmsd), median(rmsd), median(late), *after)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--duration-ps", type=int, default=10000)
    parser.add_argument("--spacing-ps", type=int, default=100)
    args = parser.parse_args()
    run = args.run.resolve()
    dest = args.output.resolve()
    dest.mkdir(exist_ok=False)
    data = []
    for case in ("no_metal", "zn"):
        for rep in (1, 2, 3):
            metrics = assess(args.gmx, run / "runs" / case / f"rep{rep}",
                             dest / f"{case}_rep{rep}", args.spacing_ps, args.duration_ps)
            data.append((case, rep, *metrics))
    with (dest / "pose_summary.tsv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")
        writer.writerow(("condition", "replicate", "frames", "final_rmsd_nm", "mean_rmsd_nm",
                         "median_rmsd_nm", "late_median_rmsd_nm", "final_c3_signed_volume_nm3",
                         "final_c4_signed_volume_nm3"))
        writer.writerows(data)


if __name__ == "__main__":
    main()
