#!/usr/bin/env python3
"""Install only the five pinned ecdata inputs, verifying every SHA-256.

Uses the Python standard library. Existing files with unexpected content are
never replaced. With --archive, the canonical archive itself is checked first.
This preparation helper is separate from the unchanged V0.2 scientific code.
"""
from argparse import ArgumentParser
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import tempfile
from urllib.request import Request, urlopen
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "25cec5ecfec8b9f016eb1631ac633194c2bed39f"
ARCHIVE_SHA256 = "1338bdf8ddb50de5d9278a14e207ee3c3f2c9ab211d929b3556d7adf92a0d4fc"
ARCHIVE_ROOT = "BSD_MW_SHAPE_STUDY_V0_2/"


def file_hash(path):
    h = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def checked_records():
    counts = json.loads((ROOT / "inputs/handoff_reference/source_counts.json").read_text())
    pins = json.loads((ROOT / "third_party/ecdata/SOURCE_FILES_SHA256.json").read_text())
    if counts["commit"] != COMMIT or pins["source_commit"] != COMMIT:
        raise ValueError("Source commit disagrees with the pinned study")
    entries = {r["path"]: r for r in pins["files"]}
    expected_blocks = ["00000-09999", "10000-19999", "20000-29999", "30000-39999", "40000-49999"]
    if [r["block"] for r in counts["blocks"]] != expected_blocks or len(entries) != 5:
        raise ValueError("Expected exactly the five study input blocks")
    records = []
    for b in counts["blocks"]:
        name = "allgens." + b["block"]
        relative = "inputs/allgens/" + name
        pin = entries[relative]
        url = f"https://raw.githubusercontent.com/JohnCremona/ecdata/{COMMIT}/allgens/{name}"
        if b["url"] != url or b["sha256"] != pin["sha256"]:
            raise ValueError(f"Inconsistent source pin for {name}")
        records.append({"path": relative, "url": url, "sha256": pin["sha256"], "bytes": pin["bytes"]})
    return records


def check_bytes(data, record):
    if len(data) != record["bytes"] or sha256(data).hexdigest() != record["sha256"]:
        raise ValueError(f"Input integrity check failed: {record['path']}")


def write_verified(data, target, record):
    check_bytes(data, record)
    target.parent.mkdir(parents=True, exist_ok=True)
    # No unverified network/archive content is installed at the final path.
    fd, name = tempfile.mkstemp(prefix=target.name + ".", dir=target.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
        if target.exists():
            raise FileExistsError(f"Refusing concurrent replacement: {target}")
        os.replace(name, target)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, help="Read the unchanged canonical V0.2 ZIP instead of downloading")
    ap.add_argument("--verify-only", action="store_true", help="Check existing inputs without downloading or writing")
    args = ap.parse_args()
    if args.archive and args.verify_only:
        ap.error("Choose --archive or --verify-only")
    records = checked_records()
    missing = []
    for r in records:
        dest = ROOT / r["path"]
        if dest.exists():
            if dest.stat().st_size != r["bytes"] or file_hash(dest) != r["sha256"]:
                raise ValueError(f"Existing input is mismatched; left untouched: {dest}")
            print("VERIFIED", r["path"])
        else:
            missing.append(r)
    if args.verify_only:
        if missing:
            raise FileNotFoundError("Missing pinned inputs: " + ", ".join(r["path"] for r in missing))
    elif args.archive:
        if file_hash(args.archive) != ARCHIVE_SHA256:
            raise ValueError("Archive SHA-256 does not match the canonical V0.2 ZIP")
        with ZipFile(args.archive) as archive:
            for r in missing:
                data = archive.read(ARCHIVE_ROOT + r["path"])
                write_verified(data, ROOT / r["path"], r)
                print("INSTALLED", r["path"])
    else:
        for r in missing:
            request = Request(r["url"], headers={"User-Agent": "MW-rank2-census-input-fetch/1"})
            with urlopen(request, timeout=60) as response:
                data = response.read(r["bytes"] + 1)
            write_verified(data, ROOT / r["path"], r)
            print("INSTALLED", r["path"])
    print(f"PASS: 5/5 input blocks verified at commit {COMMIT}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
