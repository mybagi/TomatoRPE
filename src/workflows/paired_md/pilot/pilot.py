"""Prepare and run the tomato RPE Ru5P matched pilot or a short model preflight."""

from argparse import ArgumentParser
import csv
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from tempfile import mkdtemp


PACKAGE = Path(__file__).resolve().parent
CONDITIONS = ("no_metal", "zn")
TASKS = [(condition, rep) for condition in CONDITIONS for rep in range(1, 4)]
DT_PS = 0.002
STEPS = {"nvt": 50_000, "npt": 50_000, "pilot": 5_000_000}
PREFLIGHT_PS = 500
PREFLIGHT_FRAME_PS = 10
DONORS = ((89, "NE2"), (91, "OD1"), (122, "ND1"), (233, "OD2"))
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
PROTOCOL_FILES = ("pilot.py", "reference/pocket.tsv", "summarize.py")


def pdb_atoms(path):
    atoms = []
    for line in path.read_text(encoding="ascii").splitlines():
        if line[:6].strip() in {"ATOM", "HETATM"}:
            atoms.append((
                int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                tuple(float(line[i:i+8]) for i in (30, 38, 46)),
            ))
    if len(atoms) < 100:
        raise ValueError(f"Incomplete input structure: {path}")
    return atoms


def heavy(name):
    return not name.lstrip("0123456789").upper().startswith(("H", "D"))


def fingerprint(directory):
    result = {}
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Use real files in the input snapshot: {path}")
        if path.is_file():
            with path.open("rb") as f:
                result[path.relative_to(directory).as_posix()] = hashlib_digest(f)
    return result


def hashlib_digest(file):
    checksum = sha256()
    for chunk in iter(lambda: file.read(1024 * 1024), b""):
        checksum.update(chunk)
    return checksum.hexdigest()


def inspect_condition(directory, condition, preflight=False):
    for name in ("solute.pdb", "topol.top", "model.json"):
        if not (directory / name).is_file():
            raise ValueError(f"{condition}: missing {name} in {directory}")
    metadata = json.loads((directory / "model.json").read_text(encoding="utf-8"))
    atoms = pdb_atoms(directory / "solute.pdb")
    ligand = [atom for atom in atoms if atom[1] == "MOL"]
    if len(ligand) != 23 or sum(heavy(atom[2]) for atom in ligand) != 14:
        raise ValueError(f"{condition}: expected one archived 23-atom Ru5P residue MOL")
    zinc = [i for i, atom in enumerate(atoms, 1) if atom[2].upper() in {"ZN", "ZN2"}]
    if condition == "no_metal" and zinc:
        raise ValueError("Zinc found in the no-metal reference")
    if condition == "zn":
        if len(zinc) != 1:
            raise ValueError("Expected exactly one labelled site Zn atom in the Zn model")
        reviewed = metadata.get("reviewed") is True and not metadata.get("preflight_only")
        preflight_model = preflight and metadata.get("preflight_only") is True and metadata.get("reviewed") is False
        if not (
            (reviewed or preflight_model)
            and metadata.get("substrate_unrestrained_in_production") is True
            and isinstance(metadata.get("parameter_source"), str)
            and metadata["parameter_source"].strip()
            and (directory / "model_review.md").is_file()
            and (directory / "model_review.md").stat().st_size > 0
        ):
            raise ValueError("Zn model requires parameter provenance, matching review status and an unrestrained substrate")
    return atoms, zinc


def inspect_inputs(input_root, reference_only=False, preflight=False):
    base, _ = inspect_condition(input_root / "no_metal", "no_metal", preflight)
    if reference_only:
        return
    metal, _ = inspect_condition(input_root / "zn", "zn", preflight)
    for name, atoms in (("no_metal", base), ("zn", metal)):
        coords = {(r, atom): v for r, _, atom, v in atoms if 1 <= r <= 282 and atom == "CA"}
        if len(coords) != 282:
            raise ValueError(f"{name}: residue numbering must preserve the 282 C-alpha atoms")
    lig_a = {atom: p for _, name, atom, p in base if name == "MOL" and heavy(atom)}
    lig_b = {atom: p for _, name, atom, p in metal if name == "MOL" and heavy(atom)}
    if lig_a.keys() != lig_b.keys() or any(math.dist(lig_a[k], lig_b[k]) > 0.05 for k in lig_a):
        raise ValueError("Match both Ru5P starting heavy-atom coordinates within 0.05 angstrom")
    ca_a = {r: v for r, _, a, v in base if 1 <= r <= 282 and a == "CA"}
    ca_b = {r: v for r, _, a, v in metal if 1 <= r <= 282 and a == "CA"}
    if any(math.dist(ca_a[k], ca_b[k]) > 0.05 for k in ca_a):
        raise ValueError("Match receptor C-alpha starting coordinates within 0.05 angstrom")


def prepare(inputs, work, preflight=False):
    inspect_inputs(inputs, preflight=preflight)
    work.mkdir(exist_ok=False, parents=True)
    sources = {}
    for condition in CONDITIONS:
        source = inputs / condition
        sources[condition] = fingerprint(source)
        frozen = work / "input_snapshot" / condition
        shutil.copytree(source, frozen)
        if sources[condition] != fingerprint(frozen):
            raise ValueError(f"Snapshot consistency check failed: {condition}")
    plan = {
        "preflight_only": preflight,
        "duration_ps": PREFLIGHT_PS if preflight else 10000,
        "dt_ps": DT_PS,
        "frames_every_ps": PREFLIGHT_FRAME_PS if preflight else 100,
        "input_sha256": sources,
        "protocol_sha256": {p: sha256((PACKAGE / p).read_bytes()).hexdigest() for p in PROTOCOL_FILES},
        "tasks": [
            {"condition": c, "rep": n, "gen_seed": 463101 + i, "ld_base": 584201 + 10*i}
            for i, (c, n) in enumerate(TASKS)
        ],
    }
    (work / "plan.json").write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")


def command(cwd, subcommand, *args, stdin=None, capture=False):
    cmd = [os.environ.get("GMX", "gmx"), subcommand, *map(str, args)]
    with (cwd / "commands.log").open("a", encoding="utf-8") as log:
        log.write(json.dumps(cmd) + "\n")
        log.flush()
        if capture:
            p = subprocess.run(cmd, cwd=cwd, input=stdin, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
            log.write(p.stdout)
        else:
            p = subprocess.run(cmd, cwd=cwd, input=stdin, text=True, stdout=log,
                               stderr=subprocess.STDOUT)
    p.check_returncode()
    return p.stdout if capture else ""


def write_text(path, text):
    path.write_text(text.strip() + "\n", encoding="ascii", newline="\n")


def make_index(path, groups):
    with path.open("w", encoding="ascii", newline="\n") as out:
        for name, atoms in groups.items():
            out.write(f"[ {name} ]\n")
            atoms = list(atoms)
            for i in range(0, len(atoms), 15):
                out.write(" ".join(map(str, atoms[i:i+15])) + "\n")


def gro_atoms(path):
    lines = path.read_text(encoding="ascii").splitlines()
    n = int(lines[1])
    if len(lines) != n + 3:
        raise ValueError(f"Bad atom count in {path}")
    return [(line[5:10].strip(), line[10:15].strip()) for line in lines[2:-1]]


CUTOFFS = """
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
"""


def solvate(run_dir, source, condition, preflight=False):
    setup = run_dir / "setup"
    if setup.is_dir():
        return setup
    temp = Path(mkdtemp(prefix="build_", dir=run_dir))
    shutil.copytree(source, temp, dirs_exist_ok=True)
    atoms, zinc = inspect_condition(temp, condition, preflight)
    command(temp, "--version")
    command(temp, "editconf", "-f", "solute.pdb", "-o", "box.gro", "-c",
            "-d", "1.0", "-bt", "dodecahedron")
    command(temp, "solvate", "-cp", "box.gro", "-cs", "spc216.gro", "-o",
            "solv.gro", "-p", "topol.top")
    waters = [i for i, (r, _) in enumerate(gro_atoms(temp / "solv.gro"), 1) if r == "SOL"]
    if not waters:
        raise ValueError("No SOL water atoms found for genion")
    # genion uses the selected group name to update [ molecules ] in the topology.
    make_index(temp / "water.ndx", {"SOL": waters})
    # The ion-selection TPR precedes neutralization and never integrates forces.
    ion_cutoffs = CUTOFFS.replace("coulombtype = PME", "coulombtype = Cut-off")
    write_text(temp / "ions.mdp", "integrator = steep\nnsteps = 0\n" + ion_cutoffs)
    command(temp, "grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p",
            "topol.top", "-o", "ions.tpr")
    command(temp, "genion", "-s", "ions.tpr", "-n", "water.ndx", "-o", "start.gro",
            "-p", "topol.top", "-pname", "NA", "-nname", "CL", "-conc", "0.15",
            "-neutral", "-seed", "258147", stdin="0\n")
    solvated = gro_atoms(temp / "start.gro")
    if solvated[:len(atoms)] != [(resname, name) for _, resname, name, _ in atoms]:
        raise ValueError("Solute atom order changed during setup")
    with (PACKAGE / "reference" / "pocket.tsv").open(encoding="utf-8") as f:
        pocket_resids = {int(row["resid"]) for row in csv.DictReader(f, delimiter="\t")}
    pocket_atoms = []
    ligand_atoms = []
    seen = set()
    donors = {}
    for i, (resid, resname, name, _) in enumerate(atoms, 1):
        if resname == "MOL" and heavy(name):
            ligand_atoms.append(i)
        elif resid in pocket_resids and i not in zinc and heavy(name):
            pocket_atoms.append(i)
            seen.add(resid)
        if (resid, name) in DONORS:
            donors[f"{resid}_{name}"] = i
    if seen != pocket_resids or len(seen) != 31 or len(ligand_atoms) != 14:
        raise ValueError("The fixed pocket and Ru5P atom selections do not reproduce")
    make_index(temp / "index.ndx", {
        "System": range(1, len(solvated)+1),
        "Solute": range(1, len(atoms)+1),
        "Bath": range(len(atoms)+1, len(solvated)+1),
        "PocketHeavy": pocket_atoms,
        "Ru5PHeavy": ligand_atoms,
    })
    if zinc:
        if len(donors) != len(DONORS):
            raise ValueError("Missing a predicted coordinating atom in the Zn input")
        (temp / "coordination.json").write_text(json.dumps({"zinc": zinc[0], "donors": donors}), encoding="utf-8")
    temp.rename(setup)
    return setup


def target_steps(name, preflight):
    if preflight and name == "pilot":
        return round(PREFLIGHT_PS / DT_PS)
    return STEPS[name]


def stage_mdp(name, task, preflight=False):
    if name == "em":
        return "integrator = steep\nnsteps = 50000\nemtol = 1000\nemstep = 0.01\n" + CUTOFFS
    nvt = name == "nvt"
    pressure = """pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 2.0
ref_p = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
"""
    output_steps = round(PREFLIGHT_FRAME_PS / DT_PS) if preflight else 50_000
    common = f"""integrator = md
dt = {DT_PS}
nsteps = {target_steps(name, preflight)}
tinit = 0
init-step = 0
nstxout-compressed = {output_steps}
nstenergy = 5000
nstlog = 5000
constraint-algorithm = LINCS
constraints = h-bonds
tcoupl = V-rescale
tc-grps = Solute Bath
tau_t = 0.1 0.1
ref_t = 300 300
DispCorr = EnerPres
continuation = {'no' if nvt else 'yes'}
gen-vel = {'yes' if nvt else 'no'}
ld-seed = {task['ld_base'] + ('nvt', 'npt', 'pilot').index(name)}
"""
    restraint = "define = -DPOSRES\n" if name != "pilot" else ""
    velocities = f"gen-temp = 300\ngen-seed = {task['gen_seed']}\npcoupl = no\n" if nvt else pressure
    return common + CUTOFFS + restraint + velocities


def validate_checkpoint(directory, name, preflight=False):
    dump = command(directory, "dump", "-cp", f"{name}.cpt", capture=True)
    steps = re.findall(r"(?im)^\s*step\s*=\s*(\d+)\s*$", dump)
    if not steps or int(steps[0]) < target_steps(name, preflight):
        raise ValueError(f"{name} did not reach its target step. Rerun the same task to resume.")


def run_stage(run_dir, setup, stage, previous, task, threads, preflight=False):
    stage_dir = run_dir / stage
    stage_dir.mkdir(exist_ok=True)
    source_coordinates = setup / "start.gro" if previous is None else run_dir / previous / f"{previous}.gro"
    mdp = stage_dir / "stage.mdp"
    write_text(mdp, stage_mdp(stage, task, preflight))
    top = setup / "topol.top"
    tpr = stage_dir / f"{stage}.tpr"
    if not tpr.exists():
        checkpoint_input = [] if previous in {None, "em"} else [
            "-t", run_dir / previous / f"{previous}.cpt",
        ]
        command(setup, "grompp", "-f", mdp, "-c", source_coordinates, "-r",
                setup / "start.gro", "-p", top, "-n", setup / "index.ndx",
                "-po", stage_dir / "mdout.mdp", "-pp", stage_dir / "processed.top",
                "-o", stage_dir / "pending.tpr", *checkpoint_input)
        (stage_dir / "pending.tpr").replace(tpr)

    finished = stage_dir / "finished"
    if finished.exists():
        return
    checkpoint = stage_dir / f"{stage}.cpt"
    if checkpoint.exists():
        validate = ["-cpi", checkpoint.name, "-append"]
    else:
        validate = []
        if (stage_dir / f"{stage}.log").exists():
            raise ValueError(f"{stage} already has a log but no checkpoint. Preserve that partial run.")
    command(stage_dir, "mdrun", "-s", tpr.name, "-deffnm", stage, "-ntmpi", 1,
            "-ntomp", threads, "-nb", "cpu" if stage == "em" else "gpu",
            "-cpt", "10", *validate)
    if not (stage_dir / f"{stage}.gro").is_file():
        raise ValueError(f"{stage} is missing its final coordinates")
    if stage == "em":
        text = (stage_dir / "em.log").read_text(encoding="utf-8", errors="replace")
        match = re.search(r"(?im)^\s*Maximum force\s*=\s*(" + NUMBER + ")", text)
        if not match or float(match[1]) >= 1000:
            raise ValueError("EM maximum force failed the 1000 kJ/mol/nm check")
    else:
        validate_checkpoint(stage_dir, stage, preflight)
    write_text(finished, "finished")


def analyze(run_dir, setup, preflight=False):
    from summarize import read_distances

    stage = run_dir / "pilot"
    provisional = run_dir / "pocket.pending.xvg"
    command(run_dir, "mindist", "-s", stage / "pilot.tpr", "-f", stage / "pilot.xtc",
            "-n", setup / "index.ndx", "-od", provisional, "-on", "contacts.xvg",
            "-tu", "ps", "-xvg", "none", "-pbc", stdin="3\n4\n")
    if preflight:
        read_distances(provisional, PREFLIGHT_PS // PREFLIGHT_FRAME_PS + 1, PREFLIGHT_FRAME_PS)
    else:
        read_distances(provisional)
    coordination = setup / "coordination.json"
    if coordination.exists():
        meta = json.loads(coordination.read_text(encoding="utf-8"))
        for donor, index in meta["donors"].items():
            command(run_dir, "distance", "-s", stage / "pilot.tpr", "-f", stage / "pilot.xtc",
                    "-select", f"atomnr {meta['zinc']} {index}",
                    "-oall", f"zn_{donor}_nm.xvg", "-tu", "ps", "-xvg", "none")
    provisional.replace(run_dir / "pocket.xvg")


def run(work, task_id, dedicated_threads=None):
    import fcntl

    if dedicated_threads is None:
        if not os.environ.get("SLURM_JOB_ID"):
            raise ValueError("Run workloads from an approved Slurm allocation")
        threads = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
    else:
        if os.environ.get("SLURM_JOB_ID") or shutil.which("sbatch"):
            raise ValueError("Use the scheduler entry point on a managed Slurm cluster")
        if not re.fullmatch(r"GPU-[0-9a-fA-F-]+", os.environ.get("CUDA_VISIBLE_DEVICES", "")):
            raise ValueError("Expose exactly one allocated GPU for dedicated-instance mode")
        if not 1 <= dedicated_threads <= len(os.sched_getaffinity(0)):
            raise ValueError("Dedicated-instance thread count exceeds the available CPUs")
        threads = dedicated_threads
    plan = json.loads((work / "plan.json").read_text(encoding="utf-8"))
    if plan["protocol_sha256"] != {p: sha256((PACKAGE / p).read_bytes()).hexdigest()
                                  for p in PROTOCOL_FILES}:
        raise ValueError("The frozen simulation protocol changed")
    preflight = plan.get("preflight_only") is True
    expected_schedule = (PREFLIGHT_PS, PREFLIGHT_FRAME_PS) if preflight else (10000, 100)
    if (plan["duration_ps"], plan["frames_every_ps"]) != expected_schedule:
        raise ValueError("The frozen simulation duration or output schedule is inconsistent")
    task = plan["tasks"][task_id]
    condition = task["condition"]
    source = work / "input_snapshot" / condition
    if fingerprint(source) != plan["input_sha256"][condition]:
        raise ValueError("Frozen inputs changed")
    run_dir = work / "runs" / condition / f"rep{task['rep']}"
    run_dir.mkdir(exist_ok=True, parents=True)
    with (run_dir / ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        setup = solvate(run_dir, source, condition, preflight)
        previous = None
        for stage in ("em", "nvt", "npt", "pilot"):
            run_stage(run_dir, setup, stage, previous, task, threads, preflight)
            previous = stage
        analyze(run_dir, setup, preflight)


def main():
    parser = ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    check = sub.add_parser("check")
    check.add_argument("--reference-only", action="store_true")
    provision = sub.add_parser("prepare")
    provision.add_argument("run_directory", type=Path)
    for subparser in (check, provision):
        subparser.add_argument("--inputs", type=Path, default=PACKAGE / "conditions")
        subparser.add_argument("--preflight", action="store_true",
                               help="Use a candidate model for a 500 ps validation trajectory")
    execute = sub.add_parser("run")
    execute.add_argument("run_directory", type=Path)
    execute.add_argument("task", type=int, choices=range(6))
    execute.add_argument("--dedicated-threads", type=int, help="Use one assigned GPU on an already provisioned dedicated instance")
    args = parser.parse_args()
    try:
        if args.action == "check":
            inspect_inputs(args.inputs, args.reference_only, args.preflight)
            print("Input structure/metadata check passed.")
        elif args.action == "prepare":
            prepare(args.inputs, args.run_directory, args.preflight)
            print(f"Six tasks frozen under {args.run_directory.resolve()}")
        else:
            run(args.run_directory.resolve(), args.task, args.dedicated_threads)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"pilot: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
