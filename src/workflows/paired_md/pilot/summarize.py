"""Summarize the complete pocket-distance series for six matched trajectories."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import fmean, median
import sys


CONDITIONS = ("no_metal", "zn")
REPLICATES = ("rep1", "rep2", "rep3")
FRAME_SPACING_PS = 100
FRAME_COUNT = 101
TIME_TOLERANCE_PS = 0.01
HEADER = (
    "condition", "replicate", "n_frames", "first_time_ps", "last_time_ps",
    "first_distance_nm", "last_distance_nm", "minimum_distance_nm",
    "mean_distance_nm", "median_distance_nm", "n_below_0p4_nm",
    "fraction_below_0p4_nm", "n_below_0p5_nm", "fraction_below_0p5_nm",
)


def read_distances(path: Path, frame_count: int = FRAME_COUNT,
                   frame_spacing_ps: int = FRAME_SPACING_PS) -> list[float]:
    """Read one GROMACS minimum-distance series in ps and nm."""
    distances = []
    with path.open(encoding="utf-8-sig") as source:
        for line_number, line in enumerate(source, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            try:
                time_ps, distance_nm = map(float, line.split())
            except ValueError as error:
                raise ValueError(f"line {line_number}: expected two numeric columns") from error
            if len(distances) >= frame_count:
                raise ValueError(f"line {line_number}: more than {frame_count} frames")
            expected_ps = len(distances) * frame_spacing_ps
            if not math.isclose(time_ps, expected_ps, rel_tol=0, abs_tol=TIME_TOLERANCE_PS):
                raise ValueError(f"line {line_number}: expected {expected_ps} ps, found {time_ps:g} ps")
            if not math.isfinite(distance_nm) or distance_nm < 0:
                raise ValueError(f"line {line_number}: invalid distance {distance_nm:g} nm")
            distances.append(distance_nm)
    if len(distances) != frame_count:
        end_ps = (frame_count - 1) * frame_spacing_ps
        raise ValueError(f"expected {frame_count} frames from 0 to {end_ps} ps, found {len(distances)}")
    return distances


def summarize_trajectory(distances: list[float], frame_spacing_ps: int = FRAME_SPACING_PS) -> tuple:
    below_4 = sum(distance < 0.4 for distance in distances)
    below_5 = sum(distance < 0.5 for distance in distances)
    return (
        len(distances), 0, (len(distances) - 1) * frame_spacing_ps,
        f"{distances[0]:.6f}", f"{distances[-1]:.6f}", f"{min(distances):.6f}",
        f"{fmean(distances):.6f}", f"{median(distances):.6f}",
        below_4, f"{below_4 / len(distances):.6f}", below_5, f"{below_5 / len(distances):.6f}",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="pilot directory that contains runs/")
    args = parser.parse_args()
    count, spacing = FRAME_COUNT, FRAME_SPACING_PS
    plan_file = args.run_dir / "plan.json"
    if plan_file.is_file():
        try:
            plan = json.loads(plan_file.read_text(encoding="utf-8"))
            if plan.get("preflight_only") is True:
                from pilot import PREFLIGHT_FRAME_PS, PREFLIGHT_PS
                spacing = PREFLIGHT_FRAME_PS
                count = PREFLIGHT_PS // spacing + 1
                if (plan["duration_ps"], plan["frames_every_ps"]) != (PREFLIGHT_PS, spacing):
                    raise ValueError("Invalid preflight output schedule")
        except (OSError, ValueError, KeyError) as error:
            print(f"{plan_file}: {error}", file=sys.stderr)
            return 1
    rows = []
    errors = []
    for condition in CONDITIONS:
        for replicate in REPLICATES:
            path = args.run_dir / "runs" / condition / replicate / "pocket.xvg"
            try:
                distances = read_distances(path, count, spacing)
            except (OSError, ValueError) as error:
                errors.append(f"{path}: {error}")
            else:
                rows.append((condition, replicate, *summarize_trajectory(distances, spacing)))
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADER)
    writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
