"""Continue six archived Ru5P trajectories on a dedicated three-GPU Linux instance."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timedelta, timezone
from hashlib import sha256
import json
import math
import os
from pathlib import Path
from queue import Empty, Queue
import re
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


CASES = (("no_metal", 1), ("zn", 2), ("zn", 3),
         ("zn", 1), ("no_metal", 2), ("no_metal", 3))
ARCHIVED_STEP = 5_000_000
TARGET_STEP = 50_000_000
DT_PS = 0.002
SAVED_FRAME_PS = 100
EXTENSION_PS = 90_000
GPU_COUNT = 3
PILOT_FILES = ("pilot.tpr", "pilot.cpt", "pilot.xtc", "pilot.edr", "pilot.log", "pilot.gro")


def now():
    return datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds")


def checksum(path):
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def verified_inputs(directory):
    manifest = directory / "SHA256SUMS.txt"
    hashes = {}
    for line in manifest.read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ", 1)
        relative = Path(name)
        path = directory / relative
        if (not re.fullmatch(r"[0-9a-f]{64}", expected) or name in hashes
                or relative.is_absolute() or ".." in relative.parts
                or not path.is_file() or path.resolve() != path):
            raise ValueError(f"Invalid manifest entry {name}")
        if checksum(path) != expected:
            raise ValueError(f"Checksum mismatch for {name}")
        hashes[name] = expected
    required = {"plan.json"}
    for condition, number in CASES:
        prefix = f"runs/{condition}/rep{number}/"
        required.update(prefix + "pilot/" + name for name in PILOT_FILES)
        required.update(prefix + name for name in ("setup/start.gro", "setup/index.ndx", "pocket.xvg"))
        if condition == "zn":
            required.add(prefix + "setup/coordination.json")
    if missing := required - hashes.keys():
        raise ValueError(f"Missing verified inputs: {sorted(missing)}")
    plan = json.loads((directory / "plan.json").read_text(encoding="utf-8"))
    cases = {(task["condition"], task["rep"]) for task in plan["tasks"]}
    if (plan["preflight_only"], plan["duration_ps"], plan["dt_ps"], plan["frames_every_ps"]) != (
        False, 10_000, DT_PS, SAVED_FRAME_PS
    ) or cases != set(CASES) or len(plan["tasks"]) != len(CASES):
        raise ValueError("Plan differs from the accepted six 10 ns simulations")
    return checksum(manifest), plan


def gmx_command(gmx, *args, cwd=None, env=None):
    result = subprocess.run([gmx, *map(str, args)], cwd=cwd, env=env,
                            capture_output=True, encoding="utf-8", errors="replace")
    if result.returncode:
        tail = "\n".join((result.stderr or result.stdout).splitlines()[-12:])
        raise RuntimeError(f"gmx {args[0]} exited {result.returncode}\n{tail}")
    return result.stdout + result.stderr


def verify_gmx(gmx):
    info = gmx_command(gmx, "--version")
    required = (
        r"^GROMACS version:\s+2025\.4(?:\s|$)", r"^Precision:\s+mixed\s*$",
        r"^MPI library:\s+thread_mpi\s*$", r"^GPU support:\s+CUDA\s*$",
    )
    if not all(re.search(pattern, info, re.MULTILINE) for pattern in required):
        raise RuntimeError("Use the mixed-precision GROMACS 2025.4 CUDA thread-MPI build")
    return info


def checkpoint_step(gmx, path):
    if not path.is_file():
        raise FileNotFoundError(f"Required checkpoint is missing: {path}")
    output = gmx_command(gmx, "dump", "-cp", path)
    steps = re.findall(r"^\s*step\s*=\s*(\d+)\s*$", output, re.MULTILINE)
    times = re.findall(r"^\s*t\s*=\s*(\S+)\s*$", output, re.MULTILINE)
    if len(steps) != 1 or len(times) != 1:
        raise ValueError(f"Could not read the checkpoint step and time: {path}")
    step, time_ps = int(steps[0]), float(times[0])
    if not math.isclose(time_ps, step * DT_PS, rel_tol=0, abs_tol=0.05):
        raise ValueError(f"Checkpoint time differs from the sealed 2 fs step: {path}")
    return step


def tpr_parameters(gmx, path, scratch):
    mdp = scratch / f"{path.stem}.inspect.mdp"
    # gmx dump also prints the full topology and coordinates to stdout.
    with (scratch / "dump.log").open("ab") as log:
        subprocess.run([gmx, "dump", "-s", str(path), "-om", str(mdp)],
                       cwd=scratch, stdout=subprocess.DEVNULL, stderr=log, check=True)
    values = {}
    for line in mdp.read_text(encoding="ascii").splitlines():
        key, equal, value = line.split(";", 1)[0].partition("=")
        if equal:
            values[key.strip().lower().replace("_", "-")] = value.strip()
    return values


def verify_tpr(values, nsteps):
    if (
        values["integrator"] != "md" or float(values["dt"]) != DT_PS
        or int(values["init-step"]) != 0 or float(values["tinit"]) != 0
        or int(values["nsteps"]) != nsteps
        or int(values["nstxout-compressed"]) != round(SAVED_FRAME_PS / DT_PS)
    ):
        raise ValueError("TPR differs from the archived production integrator, time, or output schedule")


def write_json(path, data):
    pending = path.with_suffix(".pending")
    pending.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    pending.replace(path)


def check_ready(case, source_id):
    receipt = json.loads((case / "resume_input.json").read_text(encoding="utf-8"))
    if (receipt["source_manifest_sha256"] != source_id
            or receipt["extended_tpr_sha256"] != checksum(case / "pilot" / "pilot_100ns.tpr")):
        raise ValueError(f"Prepared input receipt differs from current files: {case}")


def prepare(gmx, source, run_root, gmx_version):
    source_id, original_plan = verified_inputs(source)
    ready = run_root / "resume_plan.json"
    if ready.exists() and json.loads(ready.read_text(encoding="utf-8"))["source_manifest_sha256"] != source_id:
        raise ValueError("The working directory belongs to a different sealed input bundle")
    for condition, number in CASES:
        origin = source / "runs" / condition / f"rep{number}"
        target = run_root / "runs" / condition / f"rep{number}"
        target.parent.mkdir(exist_ok=True, parents=True)
        if target.exists():
            check_ready(target, source_id)
            continue
        with TemporaryDirectory(dir=target.parent, prefix=".prepare-") as pending:
            scratch = Path(pending)
            candidate = scratch / "case"
            (candidate / "pilot").mkdir(parents=True)
            shutil.copytree(origin / "setup", candidate / "setup")
            shutil.copy2(origin / "pocket.xvg", candidate / "baseline_pocket.xvg")
            for name in PILOT_FILES:
                shutil.copy2(origin / "pilot" / name, candidate / "pilot" / name)
            pilot = candidate / "pilot"
            original, updated = pilot / "pilot.tpr", pilot / "pilot_100ns.tpr"
            if checkpoint_step(gmx, pilot / "pilot.cpt") != ARCHIVED_STEP:
                raise ValueError(f"Expected the accepted 10 ns checkpoint for {origin}")
            before = tpr_parameters(gmx, original, scratch)
            verify_tpr(before, ARCHIVED_STEP)
            gmx_command(gmx, "convert-tpr", "-s", original, "-extend", EXTENSION_PS, "-o", updated, cwd=pilot)
            after = tpr_parameters(gmx, updated, scratch)
            verify_tpr(after, TARGET_STEP)
            if {k: v for k, v in before.items() if k != "nsteps"} != {
                k: v for k, v in after.items() if k != "nsteps"
            }:
                raise ValueError(f"TPR extension changed additional simulation parameters: {origin}")
            write_json(candidate / "resume_input.json", {
                "source_manifest_sha256": source_id,
                "extended_tpr_sha256": checksum(updated),
            })
            candidate.rename(target)
            print(f"Verified {condition}/rep{number}: 10 ns checkpoint and 100 ns TPR", flush=True)
    write_json(ready, {
        "source_manifest_sha256": source_id,
        "original_plan": original_plan,
        "target_ps": TARGET_STEP * DT_PS,
        "gpu_workers": GPU_COUNT,
        "prepared_at": now(),
    })
    (run_root / "gromacs-version.txt").write_text(gmx_version, encoding="utf-8")


def gpu_slots():
    devices = subprocess.run(["nvidia-smi", "--query-gpu=index,uuid,name,memory.total,driver_version",
                              "--format=csv,noheader,nounits"], check=True, capture_output=True, text=True).stdout
    slots = []
    for row in csv.reader(devices.splitlines(), skipinitialspace=True):
        if len(row) != 5 or not re.fullmatch(r"GPU-[0-9a-fA-F-]+", row[1]):
            raise ValueError(f"Unexpected GPU listing: {row}")
        slots.append(row[1])
    if len(slots) != GPU_COUNT or len(set(slots)) != GPU_COUNT:
        raise RuntimeError(f"Expected precisely {GPU_COUNT} GPUs in the dedicated instance, found {len(slots)}")
    active = subprocess.run(["nvidia-smi", "--query-compute-apps=gpu_uuid,pid",
                             "--format=csv,noheader"], check=True, capture_output=True, text=True).stdout
    if any(row.split(",")[0].strip() in slots for row in active.splitlines()):
        raise RuntimeError("A compute process is already active on the selected GPUs")
    return slots, devices


def continue_case(gmx, case, gpu, threads, start_step, record):
    pilot = case / "pilot"
    env = {**os.environ, "CUDA_DEVICE_ORDER": "PCI_BUS_ID", "CUDA_VISIBLE_DEVICES": gpu,
           "OMP_NUM_THREADS": str(threads), "LC_ALL": "C", "TZ": "Asia/Shanghai"}
    print(f"{now()} {case.parent.name}/{case.name} started on {gpu} from step {start_step}", flush=True)
    with record.open("wb") as log:
        result = subprocess.run([gmx, "mdrun", "-s", "pilot_100ns.tpr", "-deffnm", "pilot",
                                 "-cpi", "pilot.cpt", "-append", "-ntmpi", "1", "-ntomp", str(threads),
                                 "-nb", "gpu", "-cpt", "5"],
                                cwd=pilot, env=env, stdout=log, stderr=subprocess.STDOUT)
    end_step = checkpoint_step(gmx, pilot / "pilot.cpt")
    complete = result.returncode == 0 and end_step == TARGET_STEP
    outcome = "complete" if complete else "incomplete"
    if complete:
        write_json(case / "case_complete.json", {
            "end_step": end_step, "ended_at": now(), "gpu_uuid": gpu,
        })
    print(f"{now()} {case.parent.name}/{case.name}: {outcome}, step {end_step}, mdrun rc {result.returncode}",
          flush=True)
    return end_step, result.returncode, outcome


def run(gmx, run_root, threads):
    if os.environ.get("SLURM_JOB_ID") or shutil.which("sbatch"):
        raise RuntimeError("This launcher is for a dedicated instance. Use the scheduler on a managed cluster")
    plan = json.loads((run_root / "resume_plan.json").read_text(encoding="utf-8"))
    workload = Queue()
    complete = []
    for condition, number in CASES:
        case = run_root / "runs" / condition / f"rep{number}"
        check_ready(case, plan["source_manifest_sha256"])
        pilot = case / "pilot"
        for name in ("pilot.xtc", "pilot.edr", "pilot.log"):
            if not (pilot / name).is_file() or not (pilot / name).stat().st_size:
                raise ValueError(f"Missing or empty append input: {pilot / name}")
        step = checkpoint_step(gmx, pilot / "pilot.cpt")
        if not ARCHIVED_STEP <= step <= TARGET_STEP:
            raise ValueError(f"Unexpected continuation checkpoint step {step} in {case}")
        task = (condition, number, case, step)
        if step == TARGET_STEP:
            if not (case / "case_complete.json").is_file():
                write_json(case / "case_complete.json", {"end_step": step, "verified_at": now()})
            complete.append((condition, number, "", step, step, 0, "already_complete", now()))
        else:
            workload.put(task)
    if not workload.empty():
        if threads < 1 or len(os.sched_getaffinity(0)) < threads * min(GPU_COUNT, workload.qsize()):
            raise ValueError("Requested CPU threads exceed the instance affinity")
        slots, devices = gpu_slots()
        attempts = run_root / "launches"
        attempts.mkdir(exist_ok=True)
        launch = attempts / datetime.now().strftime("%Y%m%d-%H%M%S-%f")
        launch.mkdir()
        (launch / "devices.csv").write_text(devices, encoding="utf-8")

        def worker(gpu):
            status = []
            while True:
                try:
                    condition, number, case, step = workload.get_nowait()
                except Empty:
                    return status
                record = launch / f"{condition}-rep{number}.log"
                try:
                    end, rc, outcome = continue_case(gmx, case, gpu, threads, step, record)
                except (OSError, RuntimeError, ValueError) as error:
                    end, rc, outcome = "", 1, "error"
                    with record.open("a", encoding="utf-8") as log:
                        print(error, file=log)
                    print(f"{condition}/rep{number}: {error}", file=sys.stderr, flush=True)
                status.append((condition, number, gpu, step, end, rc, outcome, now()))

        with ThreadPoolExecutor(max_workers=GPU_COUNT) as pool:
            for outcomes in pool.map(worker, slots):
                complete.extend(outcomes)
    pending = run_root / "status.pending.tsv"
    with pending.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("condition", "replicate", "gpu_uuid", "start_step", "end_step",
                         "mdrun_exit_code", "status", "ended_at"))
        writer.writerows(sorted(complete))
    pending.replace(run_root / "status.tsv")
    return int(any(row[6] not in {"complete", "already_complete"} for row in complete))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("--input-root", type=Path, help="the unpacked, sealed inputs")
    parser.add_argument("--run-root", required=True, type=Path, help="a separate working directory")
    parser.add_argument("--gmx", default=os.environ.get("GMX", "gmx"))
    parser.add_argument("--threads", type=int, default=8, help="CPU threads per GPU in run mode")
    args = parser.parse_args()
    if args.action == "prepare" and args.input_root is None:
        parser.error("prepare requires --input-root")
    gmx = shutil.which(args.gmx)
    if not gmx:
        parser.error(f"GROMACS executable could not be found: {args.gmx}")
    gmx = str(Path(gmx).resolve())
    work = args.run_root.resolve()
    source = args.input_root.resolve() if args.input_root is not None else None
    if source is not None and (source == work or source in work.parents or work in source.parents):
        parser.error("use separate input and work directories")
    import fcntl

    try:
        version = verify_gmx(gmx)
        if args.action == "prepare":
            work.mkdir(exist_ok=True, parents=True)
        with (work / ".resume.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if args.action == "prepare":
                prepare(gmx, source, work, version)
                return 0
            return run(gmx, work, args.threads)
    except (OSError, RuntimeError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"tomato RPE {args.action} failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
