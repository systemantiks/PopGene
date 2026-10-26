#!/usr/bin/env python3
"""Confirm synthetic translation fixture expectations with Biopython 1.88.

Run: .venv/bin/python validation/check_translation_fixtures.py

This is an independent oracle for the expected fixture sequences used by the
Rust fixture tests. It does not compile, execute or compare production Rust.
Complete codons are translated with stops retained; incomplete end codons are
explicitly excluded. These fixtures require no alternative-initiator inference.
"""
from hashlib import sha256
import json
from pathlib import Path

import Bio
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "validation/fixtures/regression"
SUPPORTED_POPGENE_TABLES = {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}
SOURCES = [
    "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi",
    "https://www.ebi.ac.uk/ena/WebFeat/CDS_s.html",
]
EXPECTED = {
    "calculated-translations.gb": {
        "sequence": "ATGAAATAACCGTTATTTCAT",
        "frames": {
            "+1": "MK*PLFH", "+2": "*NNRYF", "+3": "EITVIS",
            "-1": "MK*RLFH", "-2": "*NNGYF", "-3": "EITVIS",
        },
        # label: (protein in coding direction, strand, codon_start, table)
        "cds": {
            "forward_MK_stop": ("MK*", 1, 1, 1),
            "reverse_MK_stop": ("MK*", -1, 1, 1),
            "unsupported_code_4": ("P", 1, 1, 4),
        },
        # These are the 0-based centers/letters asserted by the Rust layout test.
        "forward_display_glyphs": {
            "forward_MK_stop": [[1, "M"], [4, "K"], [7, "*"]],
            "reverse_MK_stop": [[13, "*"], [16, "K"], [19, "M"]],
        },
        "rust_fixture_test": (
            "crates/popgene-desktop/src/sequence_layout.rs::tests::"
            "uncached_genbank_cds_calculates_both_strands_and_additional_code"
        ),
    },
    "reading-frames.gb": {
        "sequence": "ATGAAATAACCGTTATTTCATACG",
        "frames": {
            "+1": "MK*PLFHT", "+2": "*NNRYFI", "+3": "EITVISY",
            "-1": "RMK*RLFH", "-2": "V*NNGYF", "-3": "YEITVIS",
        },
        "cds": {
            "forward_1": ("MK*PLFHT", 1, 1, 1),
            "forward_2": ("*NNRYFI", 1, 2, 1),
            "forward_3": ("EITVISY", 1, 3, 1),
            "reverse_2": ("V*NNGYF", -1, 2, 1),
        },
        "rust_fixture_test": (
            "crates/popgene-formats/src/snapgene.rs::tests::"
            "native_offset_fixture_matches_genbank_calculations_and_rebuilt_frames"
        ),
    },
}
EXPECTED["reading-frames-shifted.gb"] = {
    "sequence": EXPECTED["reading-frames.gb"]["sequence"],
    "frames": EXPECTED["reading-frames.gb"]["frames"],
    "cds": {
        "shift2_phase2": ("EITVISY", 1, 2, 1),
        "shift3_phase3": ("NNRYFI", 1, 3, 1),
        "reverse2_phase2": ("YEITVIS", -1, 2, 1),
        "reverse3_phase3": ("*NNGYF", -1, 3, 1),
    },
    "rust_fixture_test": (
        "crates/popgene-formats/src/snapgene.rs::tests::"
        "shifted_native_frames_are_feature_relative"
    ),
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def complete_translation(sequence, table):
    """Keep stops and exclude only a declared incomplete end codon."""
    complete_length = len(sequence) // 3 * 3
    return str(sequence[:complete_length].translate(table=table, to_stop=False)), len(sequence) % 3


def check_fixture(filename, expected):
    path = FIXTURES / filename
    record = SeqIO.read(path, "genbank")
    sequence = str(record.seq)
    require(sequence == expected["sequence"], f"{filename}: sequence changed")
    frames = {}
    for sign, strand in [("+", record.seq), ("-", record.seq.reverse_complement())]:
        for phase in range(3):
            protein, _ = complete_translation(strand[phase:], 1)
            frames[f"{sign}{phase + 1}"] = protein
    require(frames == expected["frames"], f"{filename}: six-frame expectations disagree")
    cds_results = []
    seen = set()
    for feature in record.features:
        if feature.type != "CDS":
            continue
        labels = feature.qualifiers.get("label", [])
        require(len(labels) == 1, f"{filename}: CDS needs one label")
        label = labels[0]
        require(label not in seen, f"{filename}: duplicate CDS {label}")
        seen.add(label)
        require(label in expected["cds"], f"{filename}: unexpected CDS {label}")
        require("translation" not in feature.qualifiers, f"{filename}: {label} has cached translation")
        starts = feature.qualifiers.get("codon_start", ["1"])
        tables = feature.qualifiers.get("transl_table", ["1"])
        require(len(starts) == len(tables) == 1, f"{filename}: conflicting CDS qualifiers")
        codon_start, table = int(starts[0]), int(tables[0])
        require(codon_start in (1, 2, 3), f"{filename}: invalid codon_start")
        phase = codon_start - 1
        # Biopython extracts biological order and reverse-complements the strand.
        coding = feature.extract(record.seq)[phase:]
        positions = list(feature.location)[phase:]
        protein, trailing_bases = complete_translation(coding, table)
        require(
            (protein, feature.location.strand, codon_start, table) == expected["cds"][label],
            f"{filename}: {label} protein/strand/codon_start/table expectations disagree",
        )
        codons = [positions[i:i + 3] for i in range(0, len(protein) * 3, 3)]
        glyphs = sorted([[codon[1], aa] for codon, aa in zip(codons, protein)])
        if label in expected.get("forward_display_glyphs", {}):
            require(glyphs == expected["forward_display_glyphs"][label],
                    f"{filename}: {label} codon center expectations disagree")
        supported = table in SUPPORTED_POPGENE_TABLES
        if not supported:
            require(label == "unsupported_code_4" and table == 4,
                    f"{filename}: unexpected unsupported table")
        cds_results.append({
            "label": label,
            "strand": feature.location.strand,
            "codon_start": codon_start,
            "transl_table": table,
            "protein_in_coding_direction": protein,
            "trailing_bases_excluded": trailing_bases,
            "codon_positions_0_based_in_coding_direction": codons,
            "oracle_glyphs_0_based_in_forward_display_order": glyphs,
            "source_translation_qualifier": None,
            "expected_popgene_behavior": (
                "display calculated protein" if supported else
                "omit translation; explicit 'table 4 is unsupported' notice; no table-1 fallback"
            ),
        })
    require(seen == set(expected["cds"]), f"{filename}: missing CDS")
    return {
        "fixture": str(path.relative_to(ROOT)),
        "sha256": sha256(path.read_bytes()).hexdigest(),
        "sequence": sequence,
        "frames_in_coding_direction": frames,
        "cds": cds_results,
        "related_rust_fixture_test": expected.get("rust_fixture_test"),
        "passed": True,
    }


def main():
    require(Bio.__version__ == "1.88", f"Expected pinned Biopython 1.88, got {Bio.__version__}")
    results = [check_fixture(filename, expected) for filename, expected in EXPECTED.items()]
    result = {
        "oracle": "Biopython 1.88",
        "scope": "Independent validation of expected synthetic fixture translations; production Rust is not executed or compared by this script.",
        "production_rust_executed": False,
        "method": "GenBank feature extraction respects strand; codon_start offsets are applied in coding direction; complete codons are translated with stops retained. Six frames use table 1. Unsupported-in-PopGene table 4 is translated only by the oracle, and its required omission notice is an expectation for separate Rust/UI tests.",
        "fixtures": results,
        "passed": True,
        "sources": SOURCES,
    }
    output = ROOT / "validation/results/calculated-translations-oracle.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"PASS: {len(results)} fixture expectation sets; Biopython {Bio.__version__}; production Rust not executed")
    print(output.relative_to(ROOT))


if __name__ == "__main__":
    main()
