"""Download and check the published GROMACS files for a system."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import tempfile
from urllib.parse import quote

from prepare_data import fetch


ASSETS = Path(__file__).resolve().parents[1] / "assets"
RECORD_URL = "https://zenodo.org/records/24425383"
RECORD_API = "https://zenodo.org/api/records/24425383"
RECORD_DOI = "10.5281/zenodo.24425383"
EPIMER_STEM = "D-Ru5P_C4-epimer"

def read_index(filename, key):
    with (ASSETS / filename).open(encoding="utf-8", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def matches_record(path, size, checksum):
    if not path.is_file() or path.stat().st_size != size:
        return False
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, lambda: hashlib.md5(usedforsecurity=False)).hexdigest()
    return digest == checksum


def published_checksums():
    with tempfile.NamedTemporaryFile(delete=False) as stream:
        path = Path(stream.name)
    try:
        fetch(RECORD_API, path)
        record = json.loads(path.read_text(encoding="utf-8"))
    finally:
        path.unlink(missing_ok=True)
    if record["doi"] != RECORD_DOI:
        raise ValueError("Zenodo returned a different record")
    return {entry["key"]: entry for entry in record["files"]}


def download_file(name, entry, output_dir):
    if Path(name).name != name or "/" in name or "\\" in name:
        raise ValueError(f"Invalid deposited filename {name!r}")
    size = int(entry["size"])
    method, checksum = entry["checksum"].split(":", maxsplit=1)
    if method != "md5":
        raise ValueError("The published checksum format changed")
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
    catalog = read_index("systems.tsv", "system")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True, choices=sorted(catalog))
    parser.add_argument("--output", required=True, type=Path, help="Directory for the three downloaded files")
    args = parser.parse_args()
    stem = EPIMER_STEM if args.system == "substrate" else args.system
    names = (f"RPE_{stem}_0ns.gro", f"RPE_{stem}_100ns.tpr", f"RPE_{stem}_100ns.xtc")
    published = published_checksums()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in names:
        print(download_file(name, published[name], args.output), flush=True)


if __name__ == "__main__":
    main()
