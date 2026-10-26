#!/usr/bin/env python3
"""Compare actual cut/rejoin annotation outputs with independent Biopython
locations, ordered extraction and translation. Covers full linear/reverse/circular
CDSs, repeated overlapping Join leaves and genuine fuzzy ends. This does not
establish public-file persistence or reference-application behavior.
"""
import json
from pathlib import Path
import subprocess

import Bio
from Bio.Restriction import EcoRI
from Bio.Seq import Seq
from Bio.SeqFeature import Location

from rust_probe import build_example

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "ATGGAATTCAAGTAA"


def span(start, end, before=False, after=False):
    return {"Range": dict(start=start, end=end, before=before, after=after)}


def location_signature(location):
    # Complement(Join(...)) and Join(Complement(...)) can encode the same
    # biological occurrence order. Compare ordered parts, strand and endpoint
    # uncertainty rather than requiring the same spelling of the syntax tree.
    return [(type(part.start).__name__, int(part.start),
             type(part.end).__name__, int(part.end), part.strand)
            for part in location.parts]


def main():
    assert Bio.__version__ == "1.88"
    # Biopython independently establishes the top-strand cut after base4.
    assert EcoRI.search(Seq(SOURCE)) == [5]
    assert tuple(map(str, EcoRI.catalyse(Seq(SOURCE)))) == (SOURCE[:4], SOURCE[4:])
    probe = build_example("popgene-core", "cloning_probe")
    cases, expectations = [], []

    def add(name, ast, expression, phase=None, reverse=False, circular=False):
        feature = dict(label=name, start=0, end=len(SOURCE), location=ast)
        if phase is not None:
            feature.update(kind="CDS", qualifiers=[["codon_start", str(phase)]])
        sequence = str(Seq(SOURCE).reverse_complement()) if reverse else SOURCE
        if circular:
            sequence = SOURCE[4:] + SOURCE[:4]
        cases.append(dict(method="digest_ligation", enzyme="EcoRI",
                          template=dict(sequence=SOURCE, circular=circular,
                                        annotations=[feature]),
                          fragments=[0] if circular else [1, 0] if reverse else [0, 1],
                          reverse=reverse, circular=circular))
        location = Location.fromstring(expression, length=len(sequence), circular=circular)
        extracted = str(location.extract(Seq(sequence)))
        coding = extracted[phase - 1:] if phase else ""
        expected_translation = str(Seq(coding[:len(coding) // 3 * 3]).translate()) if phase else None
        expectations.append(dict(name=name, sequence=sequence, circular=circular,
                                 location=expression, extracted=extracted, phase=phase,
                                 translation=expected_translation,
                                 ordered_positions=list(location)))

    for phase in (1, 2, 3):
        add(f"linear_phase{phase}", span(0, 15), "1..15", phase)
        add(f"reverse_phase{phase}", span(0, 15), "complement(1..15)", phase, reverse=True)
        add(f"circular_phase{phase}", span(0, 15), "join(12..15,1..11)", phase, circular=True)
    repeated = {"Join": [span(0, 8), span(0, 8), span(4, 12)]}
    add("ordered_repeated_join", repeated, "join(1..8,1..8,5..12)")
    add("reverse_repeated_join", repeated, "complement(join(4..11,8..15,8..15))", reverse=True)
    add("genuine_fuzzy_ends", span(0, 15, True, True), "<1..>15", 2)

    output = subprocess.run([str(probe)], input="".join(json.dumps(v) + "\n" for v in cases),
                            capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(output) == len(cases)
    for row, expected in zip(output, expectations):
        actual = json.loads(row)
        assert "error" not in actual, (expected["name"], actual)
        assert actual["sequence"] == expected["sequence"], (expected["name"], "sequence")
        assert actual["circular"] == expected["circular"], (expected["name"], "topology")
        if expected["circular"]:
            # Closed duplexes use the implicit complementary strand rather
            # than a separate physical-end record.
            assert actual["opposite"] is None
        else:
            assert actual["opposite"] == str(Seq(expected["sequence"]).reverse_complement())
        assert len(actual["annotations"]) == 1, (expected["name"], actual["annotations"])
        feature = actual["annotations"][0]
        actual_location = Location.fromstring(feature["location"], length=len(actual["sequence"]),
                                              circular=actual["circular"])
        expected_location = Location.fromstring(expected["location"], length=len(actual["sequence"]),
                                                circular=actual["circular"])
        assert location_signature(actual_location) == location_signature(expected_location), (expected["name"], feature)
        assert list(actual_location) == expected["ordered_positions"], (expected["name"], feature)
        assert feature["sequence"] == expected["extracted"], (expected["name"], feature)
        assert feature["phase"] == expected["phase"], (expected["name"], feature)
        assert feature["translation"] == expected["translation"], (expected["name"], feature)
        assert not any(k == "partial" for k, _ in feature["qualifiers"])

    report = dict(status="pass", biopython=Bio.__version__, production_cases=len(cases),
                  scope=__doc__.strip(), cases=expectations)
    output_dir = ROOT / "validation/results"
    output_dir.mkdir(exist_ok=True)
    (output_dir / "annotation-lineage-oracle.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "cases"}, indent=2))


if __name__ == "__main__":
    main()
