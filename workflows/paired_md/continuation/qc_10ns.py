"""Extract physical diagnostics from the six completed tomato RPE 10 ns trajectories.

The run tree is opened read-only. GROMACS scratch files go in a temporary
directory. The JSON report goes to stdout. Exit 2 means evidence is incomplete;
a complete report with a negative result exits 0.
"""

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
from statistics import fmean, median, pstdev
import subprocess
from tempfile import TemporaryDirectory


TASKS = [(condition, rep) for condition in ("no_metal", "zn") for rep in (1, 2, 3)]
END_PS = 10000
FRAME_PS = 100
STAGE_END_PS = {"nvt": 100, "npt": 100, "pilot": END_PS}
DONORS = ("89_NE2", "91_OD1", "122_ND1", "233_OD2")
ASPARTATE_ALTERNATES = {"91_OD1": "91_OD2", "233_OD2": "233_OD1"}
LIGAND_OXYGENS = {"Ru5P_O2_C2_carbonyl": "O2", "Ru5P_O3_C3_hydroxyl": "O3"}
NUMERICAL_WARNING = re.compile(
    r"LINCS\s+WARNING|SHAKE\s+(?:did not converge|failed)|constraint error|"
    r"(?:kinetic|potential|total) energy.{0,80}\b(?:nan|inf(?:inity)?)\b|"
    r"not finite|pressure scaling more than|particles communicated to PME node|"
    r"fatal error|segmentation fault",
    re.IGNORECASE,
)


def read_xvg(path):
    samples = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        t, value = map(float, line.split())
        if not (math.isfinite(t) and math.isfinite(value)):
            raise ValueError(f"{path.name}: non-finite data at {t:g} ps")
        if samples and t <= samples[-1][0]:
            raise ValueError(f"{path.name}: time is not strictly increasing at {t:g} ps")
        samples.append((t, value))
    if not samples:
        raise ValueError(f"{path.name}: empty time series")
    return samples


def read_distances(path):
    samples = read_xvg(path)
    times = range(0, END_PS + 1, FRAME_PS)
    if len(samples) != len(times) or any(
        abs(actual - expected) > 0.01 or distance < 0
        for expected, (actual, distance) in zip(times, samples)
    ):
        raise ValueError(f"{path.name}: expected {len(times)} nonnegative distances through {END_PS} ps")
    return samples


def describe(samples):
    values = [value for _, value in samples]
    half = len(values) // 2
    return {
        "n_samples": len(samples), "time_range_ps": [samples[0][0], samples[-1][0]],
        "first": values[0], "last": values[-1], "mean": fmean(values),
        "median": median(values), "sd": pstdev(values), "min": min(values), "max": max(values),
        "early_mean": fmean(values[:half]) if half else values[0],
        "late_mean": fmean(values[-half:]) if half else values[-1],
    }


def describe_distance(samples, cutoff, sustained_ps):
    gaps = []
    for time_ps, distance in samples:
        if distance <= cutoff:
            continue
        if gaps and abs(time_ps - gaps[-1][1] - FRAME_PS) < 0.02:
            gaps[-1][1] = time_ps
        else:
            gaps.append([time_ps, time_ps])
    longest = max((end - start for start, end in gaps), default=0)
    return {
        **describe(samples),
        "fraction_at_or_below_cutoff": sum(d <= cutoff for _, d in samples) / len(samples),
        "longest_sampled_gap_ps": longest,
        "sustained_gap_intervals_ps": [[a, b] for a, b in gaps if b - a + 0.01 >= sustained_ps],
        "samples_ps_nm": samples,
    }


def extract(gmx, scratch, command, *args, stdin=None, distance=False):
    output = scratch / "quantity.xvg"
    output.unlink(missing_ok=True)
    option = "-oall" if distance else "-o"
    result = subprocess.run(
        [gmx, command, *map(str, args), option, str(output), "-xvg", "none"],
        cwd=scratch, input=stdin or "", stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, errors="replace", env={**os.environ, "LC_ALL": "C", "GMX_NO_QUOTES": "1"},
        check=False,
    )
    if result.returncode or not output.is_file():
        detail = " ".join(result.stdout.splitlines()[-8:])
        raise ValueError(f"{command} failed ({result.returncode}): {detail}")
    return read_distances(output) if distance else read_xvg(output)


def stage_review(gmx, scratch, directory, stage, warnings):
    path = directory / f"{stage}.log"
    text = path.read_text(encoding="utf-8", errors="replace")
    warnings.extend(f"{stage}: {line.strip()}" for line in text.splitlines() if NUMERICAL_WARNING.search(line))
    if not (directory / "finished").is_file():
        raise ValueError(f"{stage}: workflow has no completion marker")
    if stage == "em":
        force_match = re.search(r"(?m)^\s*Maximum force\s*=\s*(\S+)", text)
        if not force_match:
            raise ValueError("em: maximum force is missing from the log")
        force = float(force_match[1])
        if not math.isfinite(force) or not 0 <= force < 1000:
            warnings.append(f"em: maximum force {force:g} failed 1000 kJ/mol/nm check")
        return {"maximum_force_kj_mol_nm": force if math.isfinite(force) else None}

    terms = {"Temperature": "K", "Potential": "kJ/mol"}
    if stage != "nvt":
        terms.update({"Pressure": "bar", "Density": "kg/m^3"})
    diagnostics = {}
    for term, unit in terms.items():
        samples = extract(gmx, scratch, "energy", "-f", directory / f"{stage}.edr", stdin=f"{term}\n0\n")
        if len(samples) < 2 or abs(samples[-1][0] - STAGE_END_PS[stage]) > 0.1:
            raise ValueError(f"{stage}: {term} energy series does not reach the expected end time")
        diagnostics[term] = {"units": unit, **describe(samples), "samples": samples}
    temperature = diagnostics["Temperature"]
    if abs(temperature["late_mean"] - 300) > 15:
        warnings.append(f"{stage}: late-half temperature mean is more than 15 K from the 300 K target")
    density = diagnostics.get("Density")
    if density and not 850 <= density["late_mean"] <= 1150:
        warnings.append(f"{stage}: late-half density mean is outside the broad 850–1150 kg/m^3 review range")
    return diagnostics


def zinc_series(gmx, scratch, run):
    coordination = json.loads((run / "setup" / "coordination.json").read_text(encoding="utf-8"))
    zinc = coordination["zinc"]
    atom_rows = (run / "setup" / "start.gro").read_text(encoding="utf-8").splitlines()
    atoms = atom_rows[2:2 + int(atom_rows[1])]
    donors = coordination["donors"]
    if set(donors) != set(DONORS) or not 1 <= zinc <= len(atoms) or atoms[zinc - 1][10:15].strip() not in {"ZN", "ZN2"}:
        raise ValueError("Invalid Zn site atom or expected protein donor list")

    targets = {}
    for donor in DONORS:
        atom_name = donor.split("_")[1]
        index = donors[donor]
        if not 1 <= index <= len(atoms) or atoms[index - 1][10:15].strip() != atom_name:
            raise ValueError(f"{donor}: donor index does not match the solvated structure")
        if donor in ASPARTATE_ALTERNATES:
            name = ASPARTATE_ALTERNATES[donor]
            targets[name] = (atoms[index - 1][:10], name.split("_")[1])
    for key, atom in LIGAND_OXYGENS.items():
        targets[key] = (None, atom)

    distances = {key: read_distances(run / f"zn_{key}_nm.xvg") for key in DONORS}
    stage = run / "pilot"
    for key, (residue_field, atom_name) in targets.items():
        indices = [i for i, line in enumerate(atoms, 1) if line[10:15].strip() == atom_name
                   and (line[:10] == residue_field if residue_field else line[5:10].strip() == "MOL")]
        if len(indices) != 1:
            raise ValueError(f"{key}: expected one atom in the solvated structure, found {len(indices)}")
        distances[key] = extract(gmx, scratch, "distance", "-s", stage / "pilot.tpr", "-f",
                                 stage / "pilot.xtc", "-select", f"atomnr {zinc} {indices[0]}",
                                 "-tu", "ps", "-pbc", distance=True)
    return distances


def loss_labels(metrics):
    return [name for name, values in metrics.items() if values["sustained_gap_intervals_ps"]]


def review(gmx, scratch, root, cutoff, sustained_ps, tasks=TASKS):
    reports = []
    for condition, rep in tasks:
        run = root / "runs" / condition / f"rep{rep}"
        report = {"condition": condition, "replicate": rep, "stages": {}, "warnings": [], "errors": []}
        reports.append(report)
        for stage in ("em", *STAGE_END_PS):
            try:
                report["stages"][stage] = stage_review(gmx, scratch, run / stage, stage, report["warnings"])
            except (OSError, ValueError) as error:
                report["errors"].append(f"{stage}: {error}")
        try:
            pocket = read_distances(run / "pocket.xvg")
            late_pocket = [(time, distance) for time, distance in pocket if time >= END_PS / 2]
            report["pocket_nm"] = {
                **describe(pocket), "fraction_below_0p4_nm": fmean(d < 0.4 for _, d in pocket),
                "fraction_below_0p5_nm": fmean(d < 0.5 for _, d in pocket),
                "late_n_samples": len(late_pocket),
                "late_fraction_below_0p4_nm": fmean(d < 0.4 for _, d in late_pocket),
                "late_fraction_below_0p5_nm": fmean(d < 0.5 for _, d in late_pocket),
            }
            if condition == "zn":
                distances = zinc_series(gmx, scratch, run)
                report["zn_distances_nm"] = {
                    name: describe_distance(series, cutoff, sustained_ps) for name, series in distances.items()
                }
                effective_donors = {name: distances[name] for name in DONORS}
                for name, alternate in ASPARTATE_ALTERNATES.items():
                    effective_donors[name] = [
                        (time, min(distance, alternate_distance))
                        for (time, distance), (_, alternate_distance) in zip(distances[name], distances[alternate])
                    ]
                report["zn_protein_site_nm"] = {
                    name: describe_distance(series, cutoff, sustained_ps)
                    for name, series in effective_donors.items()
                }
                report["persistent_protein_gaps"] = loss_labels(report["zn_protein_site_nm"])
                report["persistent_ligand_gaps"] = loss_labels({
                    name: report["zn_distances_nm"][name] for name in LIGAND_OXYGENS
                })
        except (OSError, ValueError, KeyError, IndexError, TypeError) as error:
            report["errors"].append(f"trajectory analysis: {error}")
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="directory with plan.json and the completed runs tree")
    parser.add_argument("--gmx", default=os.environ.get("GMX", "gmx"))
    parser.add_argument("--task", type=int, choices=range(len(TASKS)))
    parser.add_argument("--coordination-cutoff-nm", type=float, default=0.30)
    parser.add_argument("--sustained-loss-ps", type=float, default=100)
    args = parser.parse_args()
    root = args.run.resolve()
    if not (math.isfinite(args.coordination_cutoff_nm) and args.coordination_cutoff_nm > 0
            and math.isfinite(args.sustained_loss_ps) and FRAME_PS <= args.sustained_loss_ps <= END_PS):
        parser.error(f"supply a positive finite cutoff and a duration between {FRAME_PS} and {END_PS} ps")
    gmx = shutil.which(args.gmx)
    if gmx is None:
        parser.error(f"could not find GROMACS executable {args.gmx}")
    try:
        plan = json.loads((root / "plan.json").read_text(encoding="utf-8"))
        if not (plan.get("preflight_only") is False and plan["duration_ps"] == END_PS
                and plan["frames_every_ps"] == FRAME_PS
                and [(task["condition"], task["rep"]) for task in plan["tasks"]] == TASKS):
            raise ValueError("expected the six-run tomato RPE 10 ns pilot plan")
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))
    with TemporaryDirectory(prefix="rpe-10ns-qc-") as scratch:
        tasks = TASKS if args.task is None else (TASKS[args.task],)
        reports = review(str(Path(gmx).resolve()), Path(scratch), root,
                         args.coordination_cutoff_nm, args.sustained_loss_ps, tasks)
    incomplete = any(report["errors"] for report in reports)
    if incomplete:
        decision = "incomplete_evidence"
    elif any(report["warnings"] for report in reports):
        decision = "hold_for_numerical_or_thermodynamic_review"
    elif any(report.get("persistent_protein_gaps") for report in reports):
        decision = "hold_for_zinc_site_model_review"
    elif any(report.get("persistent_ligand_gaps") for report in reports):
        decision = "review_persistent_ligand_coordination_changes"
    else:
        decision = "no_persistent_coordination_gap_detected"
    print(json.dumps({
        "run": str(root), "decision": decision,
        "coordination_cutoff_nm": args.coordination_cutoff_nm,
        "sustained_loss_ps": args.sustained_loss_ps,
        "interpretation": [
            "Project triage uses the sampled gap duration. The starting six contacts are about 0.22–0.24 nm.",
            "Protein Asp site retention uses the closer of both carboxylate oxygens. The individual oxygens are retained separately.",
            "Ligand departures are observations for the unrestrained model. Check poses and alternative contacts before choosing production sampling.",
            "Pressure is reported descriptively because short trajectories have large fluctuations. Check final geometry, stable thermodynamics and possible water-mediated changes manually.",
            "When a protein donor leaves, inspect the equilibrated structure, protonation, neighboring waters and alternative donors. Repeated loss may warrant a different independently validated model before production.",
        ],
        "reports": reports,
    }, ensure_ascii=False, indent=2, allow_nan=False))
    return 2 if incomplete else 0


if __name__ == "__main__":
    raise SystemExit(main())
