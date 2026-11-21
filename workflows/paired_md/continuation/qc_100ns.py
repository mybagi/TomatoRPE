"""Review appended tomato RPE 100 ns trajectories without modifying the run tree."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
from statistics import fmean, median
import subprocess

from compare_poses import chiral_volumes, gro_atoms, indices
from qc_10ns import (
    ASPARTATE_ALTERNATES, DONORS, LIGAND_OXYGENS, NUMERICAL_WARNING, TASKS,
    describe, describe_distance, loss_labels, read_xvg,
)


BASELINE_PS = 10_000
END_PS = 100_000
LATE_PS = 50_000
FRAME_PS = 100
LAST_STEP = 50_000_000
PREFIX_TOLERANCE_NM = 1e-5
ZN_CUTOFF_NM = 0.30
SUBSTRATE_CUTOFF_NM = 0.25
SUSTAINED_PS = 100
ENERGY_TERMS = {"Temperature": "K", "Pressure": "bar", "Density": "kg/m^3", "Potential": "kJ/mol"}


def gmx_command(gmx, output, *args, stdin=""):
    command = [gmx, *map(str, args)]
    result = subprocess.run(command, input=stdin, text=True, errors="replace",
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=output,
                            env={**os.environ, "LC_ALL": "C", "GMX_NO_QUOTES": "1"}, check=False)
    with (output / "commands.log").open("a", encoding="utf-8") as log:
        print(json.dumps(command), file=log)
        if args[0] == "dump" and result.returncode == 0:
            print("\n".join(result.stdout.splitlines()[:24]), file=log)
        else:
            print(result.stdout, file=log)
    if result.returncode:
        raise ValueError(f"gmx {args[0]} exited {result.returncode}; see {output / 'commands.log'}")
    return result.stdout


def frame_series(path, end_ps=END_PS):
    samples = read_xvg(path)
    times = range(0, end_ps + 1, FRAME_PS)
    if len(samples) != len(times) or any(
        abs(t - expected) > 0.01 or value < 0 for (t, value), expected in zip(samples, times)
    ):
        raise ValueError(f"{path.name}: expected {len(times)} nonnegative values at 100 ps intervals through {end_ps} ps")
    return samples


def pocket_review(gmx, run, output, tpr, xtc):
    baseline = frame_series(run / "baseline_pocket.xvg", BASELINE_PS)
    path = output / "pocket_100ns.xvg"
    gmx_command(gmx, output, "mindist", "-s", tpr, "-f", xtc,
                "-n", run / "setup/index.ndx", "-od", path,
                "-on", output / "contacts_100ns.xvg", "-tu", "ps", "-xvg", "none", "-pbc", stdin="3\n4\n")
    samples = frame_series(path)
    deltas = [abs(old - new) for (_, old), (_, new) in zip(baseline, samples[:len(baseline)])]
    measures = {"prefix_frames": len(baseline), "max_prefix_delta_nm": max(deltas)}
    for label, window in (("all", samples), ("late", [row for row in samples if row[0] >= LATE_PS - 0.01])):
        measures[f"{label}_frames"] = len(window)
        measures[f"{label}_mean_nm"] = fmean(value for _, value in window)
        for cutoff, suffix in ((0.4, "4A"), (0.5, "5A")):
            count = sum(value < cutoff for _, value in window)
            measures[f"{label}_n_below_{suffix}"] = count
            measures[f"{label}_fraction_below_{suffix}"] = count / len(window)
    measures["last_pocket_nm"] = samples[-1][1]
    return measures


def numerical_review(gmx, stage, output, report):
    log = (stage / "pilot.log").read_text(encoding="utf-8", errors="replace")
    report["numerical_warnings"].extend(line.strip() for line in log.splitlines() if NUMERICAL_WARNING.search(line))
    checkpoint = gmx_command(gmx, output, "dump", "-cp", stage / "pilot.cpt")
    step = re.search(r"(?im)^\s*step\s*=\s*(\d+)\s*$", checkpoint)
    if not step or int(step[1]) != LAST_STEP:
        raise ValueError(f"checkpoint does not report step {LAST_STEP}")
    report["checkpoint_step"] = int(step[1])
    energies = {}
    for term, units in ENERGY_TERMS.items():
        path = output / f"{term.lower()}_100ns.xvg"
        gmx_command(gmx, output, "energy", "-f", stage / "pilot.edr", "-o", path,
                    "-xvg", "none", stdin=f"{term}\n0\n")
        samples = read_xvg(path)
        if (len(samples) < 2 or abs(samples[0][0]) > 0.1 or abs(samples[-1][0] - END_PS) > 0.1
                or any(b[0] - a[0] > 10.1 for a, b in zip(samples, samples[1:]))):
            raise ValueError(f"{term} energy is missing the start, finish or intermediate samples at the original 10 ps stride")
        late = [sample for sample in samples if sample[0] >= LATE_PS - 0.01]
        energies[term] = {"units": units, **describe(samples), "last_50ns": describe(late)}
    if abs(energies["Temperature"]["last_50ns"]["mean"] - 300) > 15:
        report["numerical_warnings"].append("Mean temperature at 50–100 ns differs by more than 15 K from 300 K")
    if not 850 <= energies["Density"]["last_50ns"]["mean"] <= 1150:
        report["numerical_warnings"].append("Mean density at 50–100 ns is outside 850–1150 kg/m^3")
    return energies


def zinc_review(gmx, run, output, tpr, xtc, atoms):
    meta = json.loads((run / "setup/coordination.json").read_text(encoding="utf-8"))
    zinc = meta["zinc"]
    if not 1 <= zinc <= len(atoms) or atoms[zinc - 1][2].upper() not in {"ZN", "ZN2"}:
        raise ValueError("recorded zinc atom is absent from the solvated structure")
    if set(meta["donors"]) != set(DONORS):
        raise ValueError("protein donor list differs from the 10 ns study")
    targets = dict(meta["donors"])
    for name, index in meta["donors"].items():
        resid, atom = name.split("_")
        if not 1 <= index <= len(atoms) or (atoms[index - 1][0], atoms[index - 1][2]) != (int(resid), atom):
            raise ValueError(f"protein donor index does not identify {name}")
    alternate_targets = {name: (int(name.split("_")[0]), "ASP", name.split("_")[1])
                         for name in ASPARTATE_ALTERNATES.values()}
    alternate_targets.update({name: (None, "MOL", atom) for name, atom in LIGAND_OXYGENS.items()})
    for name, (resid, residue, atom) in alternate_targets.items():
        matches = [index for index, (r, s, a, _) in enumerate(atoms, 1)
                   if s == residue and a == atom and (resid is None or r == resid)]
        if len(matches) != 1:
            raise ValueError(f"expected one solvated atom for {name}, found {len(matches)}")
        targets[name] = matches[0]

    distances = {}
    for name, target in targets.items():
        path = output / f"zn_{name}_nm.xvg"
        gmx_command(gmx, output, "distance", "-s", tpr, "-f", xtc, "-select", f"atomnr {zinc} {target}",
                    "-oall", path, "-tu", "ps", "-xvg", "none", "-pbc")
        distances[name] = frame_series(path)
    effective = {name: distances[name] for name in DONORS}
    for name, alternate in ASPARTATE_ALTERNATES.items():
        effective[name] = [(t, min(a, b)) for (t, a), (_, b) in zip(distances[name], distances[alternate])]

    def summaries(series):
        metrics = {}
        for key, values in series.items():
            summary = describe_distance(values, ZN_CUTOFF_NM, SUSTAINED_PS)
            del summary["samples_ps_nm"]
            summary["last_50ns"] = describe([row for row in values if row[0] >= LATE_PS - 0.01])
            metrics[key] = summary
        return metrics

    protein = summaries(effective)
    individual = summaries(distances)
    ligand = [individual[name] for name in LIGAND_OXYGENS]
    simultaneous = [(t, max(a, b) <= SUBSTRATE_CUTOFF_NM)
                    for (t, a), (_, b) in zip(*(distances[name] for name in LIGAND_OXYGENS))]
    late = [paired for time, paired in simultaneous if time >= LATE_PS - 0.01]
    return {
        "protein_donors_nm": protein, "individual_distances_nm": individual,
        "persistent_protein_gaps": loss_labels(protein),
        "persistent_substrate_gaps": [name for name, values in zip(LIGAND_OXYGENS, ligand)
                                      if values["sustained_gap_intervals_ps"]],
        "all_both_oxygen_frames": sum(paired for _, paired in simultaneous),
        "late_both_oxygen_frames": sum(late),
    }


def pose_review(gmx, output, tpr, xtc, atoms):
    initial_volumes = chiral_volumes(atoms)
    core = [i for i, (r, _, atom, _) in enumerate(atoms, 1) if 48 <= r <= 279 and atom == "CA"]
    ligand = [i for i, (_, res, a, _) in enumerate(atoms, 1) if res == "MOL" and not a.lstrip("0123456789").startswith("H")]
    solute = [i for i, (_, res, _, _) in enumerate(atoms, 1) if res not in {"SOL", "NA", "CL"}]
    if len(core) != 232 or len(ligand) != 14 or solute != list(range(1, len(solute) + 1)):
        raise ValueError("the original core, ligand or contiguous solute selection does not reproduce")
    imaging = output / "imaging.ndx"
    indices(imaging, {"CoreCa": core, "Solute": solute})
    imaged = output / "solute_centered.xtc"
    initial = output / "initial_centered.gro"
    final = output / "final_centered.gro"
    command = ("trjconv", "-s", tpr, "-f", xtc, "-n", imaging, "-pbc", "mol", "-center", "-ur", "compact")
    gmx_command(gmx, output, *command, "-o", imaged, stdin="0\n1\n")
    for time, target in ((0, initial), (END_PS, final)):
        gmx_command(gmx, output, *command, "-dump", time, "-o", target, stdin="0\n1\n")
    final_volumes = chiral_volumes(gro_atoms(final))
    if any(not math.isfinite(value) or value >= 0 for value in (*initial_volumes, *final_volumes)):
        raise ValueError("Ru5P C3/C4 geometry differs from the original native stereochemistry check")
    group_file = output / "rmsd.ndx"
    indices(group_file, {"CoreCa": core, "Ru5PHeavy": ligand})
    rmsd = output / "ligand_core_aligned_rmsd_nm.xvg"
    gmx_command(gmx, output, "rms", "-s", initial, "-f", imaged, "-n", group_file,
                "-o", rmsd, "-tu", "ps", "-xvg", "none", "-nomw", stdin="0\n1\n")
    samples = frame_series(rmsd)
    if samples[0][1] > 0.005:
        raise ValueError("imaged production origin and RMSD reference disagree")
    return {
        "all_median_rmsd_nm": median(value for _, value in samples),
        "late_median_rmsd_nm": median(value for time, value in samples if time >= LATE_PS - 0.01),
        "last_rmsd_nm": samples[-1][1],
        "initial_C3_C4_signed_volumes_nm3": initial_volumes,
        "final_C3_C4_signed_volumes_nm3": final_volumes,
    }


def review_run(gmx, root, destination, condition, rep):
    run = root / "runs" / condition / f"rep{rep}"
    output = destination / f"{condition}_rep{rep}"
    output.mkdir()
    stage = run / "pilot"
    tpr, xtc = stage / "pilot_100ns.tpr", stage / "pilot.xtc"
    report = {"condition": condition, "replicate": rep, "errors": [], "numerical_warnings": []}
    checks = {
        "pocket": lambda: pocket_review(gmx, run, output, tpr, xtc),
        "energy": lambda: numerical_review(gmx, stage, output, report),
        "pose": lambda: pose_review(gmx, output, tpr, xtc, gro_atoms(run / "setup/start.gro")),
    }
    if condition == "zn":
        checks["zinc"] = lambda: zinc_review(gmx, run, output, tpr, xtc, gro_atoms(run / "setup/start.gro"))
    for check, function in checks.items():
        try:
            report[check] = function()
        except (OSError, ValueError, AssertionError, KeyError, IndexError, TypeError, StopIteration) as error:
            report["errors"].append(f"{check}: {error}")
    if "pocket" in report and report["pocket"]["max_prefix_delta_nm"] > PREFIX_TOLERANCE_NM:
        report["errors"].append(f"pocket: {BASELINE_PS} ps prefix differs from the archived trace")
    return report


def save_reports(destination, reports):
    fields = ["condition", "replicate", "status", "prefix_frames", "max_prefix_delta_nm"]
    for label in ("all", "late"):
        fields.extend(f"{label}_{name}" for name in (
            "frames", "mean_nm", "n_below_4A", "fraction_below_4A", "n_below_5A", "fraction_below_5A"))
    fields += ["last_pocket_nm", "all_median_rmsd_nm", "late_median_rmsd_nm", "last_rmsd_nm",
               "all_both_oxygen_frames", "late_both_oxygen_frames", "persistent_protein_gaps",
               "persistent_substrate_gaps"]
    with (destination / "summary_100ns.tsv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for report in reports:
            zinc = report.get("zinc", {})
            status = "incomplete" if report["errors"] else "numerical_review" if report["numerical_warnings"] else "complete"
            row = {"condition": report["condition"], "replicate": report["replicate"], "status": status,
                   **report.get("pocket", {}), **report.get("pose", {}), **zinc}
            for field in ("persistent_protein_gaps", "persistent_substrate_gaps"):
                row[field] = ",".join(zinc.get(field, []))
            writer.writerow({key: f"{value:.10g}" if isinstance(value, float) else value for key, value in row.items()})
    (destination / "physical_review_100ns.json").write_text(json.dumps({
        "duration_ps": END_PS, "late_window_ps": [LATE_PS, END_PS], "trajectory_stride_ps": FRAME_PS,
        "pocket_prefix_tolerance_nm": PREFIX_TOLERANCE_NM, "coordination_review_cutoff_nm": ZN_CUTOFF_NM,
        "both_substrate_oxygens_cutoff_nm": SUBSTRATE_CUTOFF_NM, "coordination_sampled_gap_ps": SUSTAINED_PS,
        "interpretation": [
            "Pocket proximity and Ru5P pose differences are distinct. Replicates are independent, frames are trajectory observations.",
            "Asp retention uses the nearer carboxylate oxygen. The two oxygen traces remain accessible separately.",
            "Coordination threshold excursions are observations for this Zn model and are not run failures.",
        ], "reports": reports,
    }, allow_nan=False, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--gmx", required=True)
    args = parser.parse_args()
    root, destination = args.run_root.resolve(), args.output_dir.resolve()
    if not (root / "runs").is_dir():
        parser.error(f"run directory {root / 'runs'} does not exist")
    if destination.is_relative_to(root / "runs") or destination == root:
        parser.error("choose a new output directory outside the source trajectories")
    gmx = shutil.which(args.gmx)
    if gmx is None:
        parser.error(f"GROMACS executable is unavailable at {args.gmx}")
    destination.mkdir(parents=True, exist_ok=False)
    reports = [review_run(str(Path(gmx).resolve()), root, destination, condition, rep) for condition, rep in TASKS]
    save_reports(destination, reports)
    for report in reports:
        for note in report["errors"] + report["numerical_warnings"]:
            print(f"{report['condition']}/rep{report['replicate']}: {note}")
    print(destination / "summary_100ns.tsv")
    return 2 if any(report["errors"] or report["numerical_warnings"] for report in reports) else 0


if __name__ == "__main__":
    raise SystemExit(main())
