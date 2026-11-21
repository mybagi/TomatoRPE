"""Download and check the published GROMACS files for a system."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.parse import quote


ASSETS = Path(__file__).resolve().parents[1] / "assets"
RECORD_URL = "https://zenodo.org/records/24425383"
FILE_COLUMNS = ("initial_structure", "run_input", "trajectory")
CURL = (
    "curl", "--fail", "--silent", "--show-error", "--location", "--proto", "=https",
    "--proto-redir", "=https", "--connect-timeout", "15", "--retry", "2",
)


def read_index(filename, key):
    with (ASSETS / filename).open(encoding="utf-8", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def matches_record(path, size, checksum):
    if not path.is_file() or path.stat().st_size != size:
        return False
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, lambda: hashlib.md5(usedforsecurity=False)).hexdigest()
    return digest == checksum


def fetch(url, destination):
    if os.name == "nt":
        shell = shutil.which("pwsh") or shutil.which("powershell")
        if shell is None:
            raise RuntimeError("PowerShell is required for automatic downloads on Windows")
        env = {
            **os.environ,
            "TOMATORPE_DOWNLOAD_URL": url,
            "TOMATORPE_DOWNLOAD_PATH": str(destination),
        }
        command = [
            shell, "-NoProfile", "-NonInteractive", "-Command",
            "$ProgressPreference='SilentlyContinue'; "
            "Invoke-WebRequest -UseBasicParsing -ErrorAction Stop "
            "-Uri $env:TOMATORPE_DOWNLOAD_URL -OutFile $env:TOMATORPE_DOWNLOAD_PATH",
        ]
        subprocess.run(command, env=env, check=True)
    else:
        subprocess.run(
            [*CURL, "--speed-time", "45", "--speed-limit", "512", "--output", str(destination), url],
            check=True,
        )


def download_file(name, entry, output_dir):
    if Path(name).name != name or "/" in name or "\\" in name:
        raise ValueError(f"Invalid deposited filename {name!r}")
    size = int(entry["bytes"])
    checksum = entry["md5"]
    destination = output_dir / name
    if destination.exists():
        if matches_record(destination, size, checksum):
            return destination
        raise ValueError(f"The existing file does not match Zenodo. Choose another output directory for {name}")

    url = f"{RECORD_URL}/files/{quote(name, safe='')}"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", dir=output_dir, delete=False) as stream:
            temporary = Path(stream.name)
        fetch(url, temporary)
        if not matches_record(temporary, size, checksum):
            raise ValueError(f"The downloaded file does not match the published checksum for {name}")
        temporary.replace(destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return destination


def main():
    catalog = read_index("trajectory_systems.tsv", "system")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True, choices=sorted(catalog))
    parser.add_argument("--output", required=True, type=Path, help="Directory for the three downloaded files")
    args = parser.parse_args()
    names = [catalog[args.system][column] for column in FILE_COLUMNS]
    published = read_index("trajectory_files.tsv", "file")
    args.output.mkdir(parents=True, exist_ok=True)
    for name in names:
        print(download_file(name, published[name], args.output), flush=True)


if __name__ == "__main__":
    main()
