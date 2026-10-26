#!/usr/bin/env python3
"""Check actual production figure exports with Poppler and standard PNG parsing.

Generate the inputs first:
  POPGENE_FIGURE_QA_DIR=$PWD/validation/results/figure-qa \
    scripts/cargo-local test --locked -p popgene-desktop figure

This validates dimensions, PNG resolution/CRC, the typography fixture's actual
PDF labels, and every paginated sequence row exactly once in PDF text extraction.
Native dialog/visual acceptance is separate.
"""
from pathlib import Path
import json
import re
import struct
import subprocess
import xml.etree.ElementTree as ET
import zlib

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation/results/figure-qa"


def run(*args):
    return subprocess.check_output(args, text=True, encoding="utf-8")


def check_typography_text(path):
    # These are the visible labels in figure.rs's production export fixture.
    # Checking vector paths or PDF dimensions alone also accepts an export in
    # which font lookup silently removed every text node.
    expected = "α DNA 5′ → 3′"
    extracted = run("pdftotext", "-enc", "UTF-8", "-layout", str(path), "-")
    normalized = " ".join(extracted.split())
    assert normalized == expected, (
        f"Typography labels missing or changed in {path}: "
        f"expected {expected!r}, extracted {normalized!r}"
    )
    return normalized


def main():
    info = run("pdfinfo", str(OUT / "typography.pdf"))
    size = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", info)
    assert size and tuple(map(float, size.groups())) == (300., 225.), info
    typography = check_typography_text(OUT / "typography.pdf")
    data = (OUT / "typography.png").read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    assert struct.unpack(">II", data[16:24]) == (800, 600)
    offset = 8
    resolution = None
    while offset < len(data):
        length = struct.unpack(">I", data[offset:offset + 4])[0]
        kind = data[offset + 4:offset + 8]
        payload = data[offset + 8:offset + 8 + length]
        crc = struct.unpack(">I", data[offset + 8 + length:offset + 12 + length])[0]
        assert zlib.crc32(kind + payload) == crc
        if kind == b"pHYs":
            resolution = struct.unpack(">IIB", payload)
        offset += length + 12
    assert resolution == (7559, 7559, 1), resolution
    root = ET.parse(OUT / "paginated-sequence.svg").getroot()
    expected = [int(node.attrib["data-start"]) for node in root.iter()
                if "data-popgene-row-start" in node.attrib]
    text = run("pdftotext", "-enc", "UTF-8", "-layout", str(OUT / "paginated-sequence.pdf"), "-")
    starts = [int(line.split()[0]) for line in text.splitlines()
              if re.match(r"^\d+\s+\d+", line)]
    assert starts == expected, (expected, starts)
    pages = int(re.search(r"Pages:\s+(\d+)", run("pdfinfo", str(OUT / "paginated-sequence.pdf"))).group(1))
    assert pages > 1 and text.count("Page ") == pages
    result = {"status": "pass", "pdf_dimensions_points": [300, 225],
              "typography_pdf_text": typography,
              "png_pixels": [800, 600], "png_pixels_per_metre": 7559,
              "sequence_pages": pages, "rows_exactly_once_in_order": len(expected),
              "scope": "Production vector/raster outputs checked with Poppler; native UI and visual review are separate."}
    (ROOT / "validation/results/figure-export-check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
