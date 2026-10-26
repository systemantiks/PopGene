#!/usr/bin/env python3
"""Generate minimal native regression files from reviewed synthetic records.

Only header, sequence, feature, primer and notes blocks are accepted. The source
contains no opaque history, application preferences or vendor enzyme catalogue.
"""
import argparse
import json
from pathlib import Path
import struct
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent / "fixtures"


def generate(check=False):
    source = json.loads((ROOT / "native-records.json").read_text())
    for name, records in source["files"].items():
        target = ROOT / name
        if target.resolve().parent not in (ROOT.resolve(), (ROOT / "regression").resolve()):
            raise ValueError(f"Unexpected fixture path: {name}")
        output = bytearray()
        for record in records:
            kind = record["kind"]
            if kind == 9:
                body = bytes.fromhex(record["hex"])
            elif kind == 0:
                body = bytes([record["flags"]]) + record["sequence"].encode("ascii")
            elif kind in (5, 6, 10):
                body = record["xml"].encode("utf-8")
                xml = ET.fromstring(body)
                if any(node.tag in ("UUID", "Created", "LastModified") for node in xml.iter()):
                    raise ValueError(f"Unnecessary source metadata in {name}")
            else:
                raise ValueError(f"Opaque native block {kind} is not allowed")
            output.extend(struct.pack(">BI", kind, len(body)))
            output.extend(body)
        if check:
            if target.read_bytes() != output:
                raise SystemExit(f"Regenerate fixture: {name}")
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(output)
    print(f"{'Checked' if check else 'Generated'} {len(source['files'])} synthetic native fixtures")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    generate(parser.parse_args().check)
