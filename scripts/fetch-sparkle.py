#!/usr/bin/env python3
"""Fetch the pinned upstream Sparkle binary distribution, verifying its digest."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def fetch():
    spec = json.loads((ROOT / "scripts/sparkle.json").read_text())
    cache = ROOT / ".tools/sparkle"
    cache.mkdir(parents=True, exist_ok=True)
    archive = cache / ("Sparkle-" + spec["version"] + ".tar.xz")
    if not archive.exists():
        with tempfile.TemporaryDirectory(dir=cache) as work:
            download = Path(work) / "download.tar.xz"
            subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
                            "--proto", "=https", "--proto-redir", "=https",
                            spec["url"], "-o", str(download)], check=True)
            verify(download, spec["sha256"])
            download.replace(archive)
    verify(archive, spec["sha256"])
    destination = cache / spec["version"]
    marker = destination / ".verified-sha256"
    if not marker.exists() or marker.read_text().strip() != spec["sha256"]:
        with tempfile.TemporaryDirectory(dir=cache) as work:
            unpacked = Path(work) / "distribution"
            unpacked.mkdir()
            subprocess.run(["tar", "-xJf", str(archive), "-C", str(unpacked)], check=True)
            (unpacked / ".verified-sha256").write_text(spec["sha256"] + "\n")
            if destination.exists():
                shutil.rmtree(destination)
            unpacked.replace(destination)
    return destination


def verify(path, expected):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit("Sparkle archive checksum mismatch; refusing to use " + str(path))


if __name__ == "__main__":
    print(fetch())
