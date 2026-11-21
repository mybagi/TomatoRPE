"""Prepare the reference tables used by the summary and trajectory checks."""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile


DATA_DIR = Path(__file__).resolve().parent / ".reference-data"
ARCHIVE_URL = "https://github.com/mybagi/TomatoRPE/releases/download/reference-v1/tomatorpe-reference-data.zip"
ARCHIVE_SHA256 = "192b901bd923b54214857ee3f04206c9e3f1e962617a4fde0039e26ac0676834"
CURL = (
    "curl", "--fail", "--silent", "--show-error", "--location", "--proto", "=https",
    "--proto-redir", "=https", "--connect-timeout", "15", "--retry", "2",
)


def read_rows(path, delimiter="\t"):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def verify_files(data):
    manifest = read_rows(data / "manifest.tsv")
    for entry in manifest:
        path = (data / entry["path"]).resolve()
        if not path.is_relative_to(data.resolve()):
            raise ValueError(f"Invalid manifest path: {entry['path']}")
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        if path.stat().st_size != int(entry["bytes"]) or digest != entry["sha256"]:
            raise ValueError(f"Archive integrity failure: {entry['path']}")
    return len(manifest)


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


def prepare_reference(output, archive=None):
    if output.exists():
        verify_files(output)
        return output

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="tomatorpe-data-", dir=output.parent) as folder:
        temporary = Path(folder)
        if archive is None:
            archive = temporary / "reference.zip"
            fetch(ARCHIVE_URL, archive)
        with archive.open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != ARCHIVE_SHA256:
                raise ValueError("The reference archive does not match the published checksum")

        staged = temporary / "tables"
        with zipfile.ZipFile(archive) as archive_file:
            archive_file.extractall(staged)
        verify_files(staged)
        try:
            staged.rename(output)
        except OSError:
            if not output.is_dir():
                raise
            verify_files(output)
    return output


def reference_directory(path):
    if path is None:
        return prepare_reference(DATA_DIR)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Path to a downloaded reference archive")
    parser.add_argument("--output", type=Path, default=DATA_DIR, help="Directory for the extracted tables")
    args = parser.parse_args()
    print(prepare_reference(args.output.resolve(), args.archive))


if __name__ == "__main__":
    main()
