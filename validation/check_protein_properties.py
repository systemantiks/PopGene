#!/usr/bin/env python3
"""Cross-check production Rust protein properties with Biopython 1.88.

Run: .venv/bin/python validation/check_protein_properties.py
No SnapGene result is used as the mathematical oracle. Default pI search bounds
are widened to 0–14 in both implementations; no value is accepted without a root.
"""
from collections import Counter
import json
from pathlib import Path
import random
import subprocess

import Bio
from rust_probe import build_example
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parents[1]
STANDARD = "ACDEFGHIKLMNPQRSTVWY"


def main():
    canonical = "".join((ROOT / "validation/fixtures/regression/protein-canonical.faa")
                        .read_text().splitlines()[1:])
    rng = random.Random(823)
    cases = [canonical, "A", "AA", "GG", "PETER", "INGAR", "WYC", "WYCCC",
             "D" * 100, "K" * 100, "D" * 100_000]
    cases += ["".join(rng.choices(STANDARD, k=n)) for n in (2, 3, 10, 50, 250, 5000)]
    unsupported = ["", "AX", "AB", "AZ", "AJ", "AU", "AO", "A*", "A*A", "ACGMKWVTF*XZ"]
    binary = build_example('popgene-core', 'protein_probe')
    rows = subprocess.run([str(binary)], input="\n".join(cases + unsupported) + "\n", text=True, capture_output=True, check=True).stdout.splitlines()
    assert len(rows) == len(cases + unsupported)
    for sequence, row in zip(cases + unsupported, rows):
        values = row.split("\t")
        aa_counts = Counter(sequence)
        assert int(values[6]) == len(sequence) - aa_counts["*"]
        assert int(values[7]) == aa_counts["*"]
        for cell in values[8].split(","):
            aa, count, percent = cell.split(":")
            assert int(count) == aa_counts[aa]
            expected = 100 * aa_counts[aa] / len(sequence) if sequence else 0
            assert abs(float(percent) - expected) < 1e-9
        if sequence in unsupported:
            assert values[:6] == ["-"] * 6, sequence
            continue
        analysis = ProteinAnalysis(sequence)
        assert abs(float(values[0]) - molecular_weight(sequence, "protein")) < 1e-5
        assert abs(float(values[1]) - molecular_weight(sequence, "protein", monoisotopic=True)) < 1e-5
        pi = IsoelectricPoint(sequence)
        if pi.charge_at_pH(0) < 0 or pi.charge_at_pH(14) > 0:
            assert values[2] == "-"
        else:
            assert abs(float(values[2]) - pi.pi(min_=0, max_=14)) < 1e-4
        assert abs(float(values[3]) - pi.charge_at_pH(7)) < 1e-7
        assert tuple(map(float, values[4:6])) == analysis.molar_extinction_coefficient()
    result = {
        "biopython_version": Bio.__version__,
        "canonical_cases": len(cases),
        "unsupported_cases": len(unsupported),
        "status": "pass",
        "pI_method": "Bjellqvist pKa set, pH 0–14 with root check",
        "canonical_sequence": canonical,
        "canonical_result": dict(zip(["average_da", "monoisotopic_da", "pI", "charge_pH7",
                                      "extinction_reduced", "extinction_paired_cysteines"],
                                     map(float, rows[0].split("\t")[:6]))),
    }
    (ROOT / "validation/results").mkdir(parents=True, exist_ok=True)
    (ROOT / "validation/results/protein-properties-oracle.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
