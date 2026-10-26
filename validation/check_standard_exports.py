"""Independent Biopython checks for PopGene's public-format writers.

Generate files first:
  scripts/cargo-local run -p popgene-formats --example interoperable_exports -- /tmp/popgene-exports
Then:
  .venv/bin/python validation/check_standard_exports.py /tmp/popgene-exports
Biopython is validation-only; PopGene's runtime remains Rust.
"""
import json
import sys
import warnings
from pathlib import Path

from Bio import BiopythonParserWarning, SeqIO

warnings.simplefilter("error", BiopythonParserWarning)
root = Path(__file__).resolve().parents[1]
output = Path(sys.argv[1])
source = SeqIO.read(root / "validation/fixtures/annotated.gb", "genbank")
checks = []
for name, fmt in [
    ("circular.gb", "genbank"), ("circular.embl", "embl"),
    ("circular.dna", "snapgene"), ("circular.fasta", "fasta"),
    ("protein.gp", "genbank"), ("protein.faa", "fasta"),
    ("rna.gb", "genbank"), ("rna.embl", "embl"), ("rna.frna", "fasta"),
]:
    record = SeqIO.read(output / name, fmt)
    if name.startswith("circular"):
        assert record.seq == source.seq, name
        if fmt != "fasta":
            assert record.annotations["topology"] == "circular", name
            assert len(record.features) == len(source.features), name
            # Compare by labels because native SnapGene stores primers separately.
            actual = {f.qualifiers["label"][0]: f for f in record.features}
            for expected in source.features:
                feature = actual[expected.qualifiers["label"][0]]
                assert feature.location == expected.location, (name, feature.location)
                assert feature.extract(record.seq) == expected.extract(source.seq), name
                if fmt in ("genbank", "embl"):
                    assert feature.qualifiers == expected.qualifiers, name
        else:
            assert not record.features
    if name.startswith("protein"):
        assert str(record.seq) == "ACGMKWVTF*XZ", name
        if fmt == "genbank":
            assert record.annotations["molecule_type"] == "protein"
            assert record.features[0].qualifiers["label"] == ["region"]
            assert str(record.features[0].extract(record.seq)) == "GMKWVT"
    if name.startswith("rna"):
        assert str(record.seq) == "ACGUNRY", name
        if fmt != "fasta":
            assert "RNA" in record.annotations["molecule_type"]
    checks.append({"file": name, "reader": fmt, "length": len(record), "features": len(record.features), "status": "pass"})
trace = (output / "trace.ab1").read_bytes()
assert trace == (root / "validation/fixtures/synthetic-complete.ab1").read_bytes()
trace_record = SeqIO.read(output / "trace.ab1", "abi")
assert str(trace_record.seq) == "ACGTN"
assert trace_record.letter_annotations["phred_quality"] == [40, 35, 30, 20, 5]
checks.append({"file": "trace.ab1", "reader": "abi", "status": "pass", "byte_identical": True})
print(json.dumps({"checks": checks, "parser_warnings": 0}, indent=2))
