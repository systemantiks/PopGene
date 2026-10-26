"""Independent-checker regressions; no Rust build or saved fixture is needed.

Run: .venv/bin/python -m unittest discover -s validation -p test_physical_fragments.py
Synthetic contexts deliberately use an unverified digest; these tests cannot
stand in for production activation or ancestry-reconstruction tests.
"""

import copy
import io
import json
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import Reference, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import check_physical_fragments as checker


def fragment():
    reference = dict.fromkeys(checker.shapes.REFERENCE_KEYS, None)
    reference.update(title="Current fragment citation", extra=[])
    return {
        "ends": {
            "left": {"kind": "Native", "overhang": "", "enzymes": []},
            "right": {"kind": "FivePrime", "overhang": "AATT", "enzymes": ["EcoRI"]},
            "opposite_sequence": "AATTCCAT",
        },
        "duplex_annotations": {
            "features": [{"id": 2, "label": "cross boundary", "kind": "CDS",
                          "parts": [{"start": 0, "end": 8}], "reverse": False,
                          "location": {"Range": {"start": 0, "end": 8, "before": False, "after": True}},
                          "qualifiers": [["codon_start", "2"], ["citation", "[1]"]], "edited": True}],
            "sequence_colors": [{"range": {"start": 4, "end": 8}, "reverse": True, "color_rgb": 0x993366}],
        },
        "sequence_colors": [{"range": {"start": 0, "end": 2}, "reverse": False, "color_rgb": 0x112233}],
        "sequence_qualities": [0, None, 22, 93],
        "opposite_sequence_qualities": [None, 10, 11, 12, 13, 14, 15, 16],
        "methylation": {"dam": False, "dcm": True, "eco_ki": False},
        "references": [reference],
    }


def record(top="ATGG"):
    result = SeqRecord(Seq(top), id="synthetic", name="synthetic", description="Independent fixture")
    result.annotations = {"molecule_type": "DNA", "topology": "linear"}
    reference = Reference()
    reference.title = "Current fragment citation"
    result.annotations["references"] = [reference]
    result.features = [SeqFeature(SimpleLocation(0, 4), type="misc_feature", qualifiers={"label": ["ordinary"]})]
    return result


def serialize(record, format_name):
    stream = io.StringIO()
    SeqIO.write(record, stream, format_name)
    return stream.getvalue().encode()


def block(value, raw_edit=None):
    raw = json.dumps(value, ensure_ascii=True, separators=(",", ":"))
    raw = "".join(f"\\u{ord(c):04x}" if c.isspace() else c for c in raw)
    if raw_edit:
        raw = raw_edit(raw)
    framed = f"bytes={len(raw.encode())};sha256={checker.sha(raw.encode())};json={raw}"
    return ("[PopGene physical fragment v1 begin]\n"
            + "\n".join(framed[i:i+48] for i in range(0, len(framed), 48))
            + "\n" + checker.END)


def fixture(format_name="genbank", payload=None, record_edit=None, envelope_edit=None, comment_edit=None, raw_edit=None):
    standard = record()
    if record_edit:
        record_edit(standard)
    first = serialize(standard, format_name)
    _, _, records, _ = checker.ordinary_records(first)
    value = {"format": format_name, "records_sha256": checker.sha(records),
             "context_sha256": "0" * 64, "fragment": copy.deepcopy(payload if payload is not None else fragment())}
    if envelope_edit:
        envelope_edit(value)
    comment = block(value, raw_edit)
    standard.annotations["comment"] = comment_edit(comment) if comment_edit else comment
    return serialize(standard, format_name)


class PhysicalFragmentCheckerTests(unittest.TestCase):
    def check_rejected(self, data):
        with self.assertRaises((ValueError, AssertionError, TypeError, KeyError)):
            checker.inspect(data)

    def test_valid_typed_metadata_and_honest_scope(self):
        for format_name in ("genbank", "embl"):
            with self.subTest(format=format_name):
                report, _, parsed = checker.inspect(fixture(format_name))
                self.assertTrue(report["framing_schema_raw_records_and_geometry_verified"])
                self.assertFalse(report["canonical_context_digest_verified"])
                self.assertFalse(report["complete_location_projection_verified"])
                self.assertFalse(report["ancestry_restoration_verified"])
                self.assertFalse(report["authenticity_claimed"])
                self.assertEqual(parsed, fragment())
                self.assertEqual(report["geometry"], {"extent": 8, "top": {"start": 0, "end": 4},
                                                      "bottom": {"start": 0, "end": 8}})
                # A different well-formed context hash is deliberately not
                # claimed as independently verified by this bounded checker.
                other, _, _ = checker.inspect(fixture(format_name, envelope_edit=lambda v: v.update(context_sha256="f"*64)))
                self.assertFalse(other["canonical_context_digest_verified"])

    def test_framing_checksum_and_duplicate_json_keys_fail(self):
        edits = {
            "unsupported version": lambda text: text.replace("v1 begin", "v99 begin"),
            "missing end": lambda text: text.replace(checker.END, ""),
            "duplicate block": lambda text: text + "\n" + text,
            "reversed markers": lambda text: checker.END + text.replace(checker.END, ""),
            "bad checksum": lambda text: text.replace("sha256=", "sha256=0", 1),
            "bad length": lambda text: text.replace("bytes=", "bytes=9", 1),
            "nonascii wrapping": lambda text: text.replace(";json=", ";json=\u00a0", 1),
        }
        for format_name in ("genbank", "embl"):
            for name, edit in edits.items():
                with self.subTest(format=format_name, mutation=name):
                    self.check_rejected(fixture(format_name, comment_edit=edit))
            # Checksum is recomputed, so rejection proves duplicate-key checking.
            self.check_rejected(fixture(format_name, raw_edit=lambda s: s.replace('"methylation":', '"methylation":null,"methylation":', 1)))

    def test_raw_record_changes_are_rejected_without_semantic_rematching(self):
        for format_name in ("genbank", "embl"):
            original = fixture(format_name)
            changed = original.replace(b'/label="ordinary"', b'/label="changed"')
            self.assertNotEqual(original, changed)
            with self.assertRaisesRegex(ValueError, "raw annotation record hash"):
                checker.inspect(changed)

    def test_comment_folding_and_json_depth_match_the_stated_framing_scope(self):
        def split(text):
            return text.replace("fragment v1", "fragment\nv1").replace("fragment end]", "fragment\nend]")

        # GenBank preserves explicit comment newlines; the exact v1 marker no
        # longer exists. EMBL folds its CC words, restoring the known markers.
        self.check_rejected(fixture("genbank", comment_edit=split))
        report, _, _ = checker.inspect(fixture("embl", comment_edit=split))
        self.assertTrue(report["framing_schema_raw_records_and_geometry_verified"])
        nested = fixture(raw_edit=lambda text: '{"nested":'*97 + text + '}'*97)
        with self.assertRaisesRegex(AssertionError, "depth exceeds"):
            checker.inspect(nested)

    def test_nested_fields_scalars_and_occurrence_association_are_checked(self):
        mutations = {
            "unknown fragment": lambda f: f.update(Document={"sequence": "ATGG"}),
            "unknown ends": lambda f: f["ends"].update(ancestor_sequence="ATGG"),
            "null optional chemistry": lambda f: f["ends"].update(left_modification=None),
            "unknown nested end": lambda f: f["ends"]["left"].update(extra=True),
            "unknown feature": lambda f: f["duplex_annotations"]["features"][0].update(sequence="ATGG"),
            "wrong occurrence": lambda f: f["duplex_annotations"]["features"][0].update(id=1),
            "unknown nested location": lambda f: f["duplex_annotations"]["features"][0]["location"]["Range"].update(other=1),
            "bool coordinate": lambda f: f["duplex_annotations"]["features"][0]["parts"][0].update(start=False),
            "unknown reference": lambda f: f["references"][0].update(ancestor=True),
            "reference count": lambda f: f.update(references=[]),
            "citation range": lambda f: f["duplex_annotations"]["features"][0]["qualifiers"].append(["citation", "2"]),
            "unknown methylation": lambda f: f["methylation"].update(cpg=True),
            "integer boolean": lambda f: f["methylation"].update(dam=0),
            "unknown color": lambda f: f["sequence_colors"][0].update(alpha=255),
        }
        for name, mutate in mutations.items():
            with self.subTest(mutation=name):
                value = fragment()
                mutate(value)
                self.check_rejected(fixture(payload=value))

    def test_geometry_qualities_and_full_axis_bounds_are_independent_checks(self):
        mutations = {
            "opposite mismatch": lambda f: f["ends"].update(opposite_sequence="AATTCGAT"),
            "wrong overhang bases": lambda f: f["ends"]["right"].update(overhang="AAAC"),
            "wrong polarity": lambda f: f["ends"]["right"].update(kind="ThreePrime"),
            "blunt overhang": lambda f: f["ends"]["right"].update(kind="Blunt"),
            "unpaired": lambda f: f["ends"]["left"].update(kind="FivePrime", overhang="ATGG"),
            "short qualities": lambda f: f.update(sequence_qualities=[0]),
            "bool quality": lambda f: f["opposite_sequence_qualities"].__setitem__(0, False),
            "high quality": lambda f: f["sequence_qualities"].__setitem__(0, 94),
            "out of bounds full axis": lambda f: f["duplex_annotations"]["features"][0]["parts"][0].update(end=9),
            "out of bounds location": lambda f: f["duplex_annotations"]["features"][0]["location"]["Range"].update(end=9),
            "oversized gap": lambda f: f["duplex_annotations"]["features"][0].update(location={"Gap": {"length": 10_000_001, "unknown100": False}}),
            "out of bounds color": lambda f: f["duplex_annotations"]["sequence_colors"][0]["range"].update(end=9),
            "rgb overflow": lambda f: f["sequence_colors"][0].update(color_rgb=0x1000000),
        }
        for name, mutate in mutations.items():
            with self.subTest(mutation=name):
                value = fragment()
                mutate(value)
                self.check_rejected(fixture(payload=value))
        self.check_rejected(fixture().replace(b"linear", b"circular", 1))

    def test_all_four_cohesive_terminal_polarities(self):
        axis = "TTAACCGGAATTCTGGAACCGGTT"
        for left in ("FivePrime", "ThreePrime"):
            for right in ("FivePrime", "ThreePrime"):
                with self.subTest(left=left, right=right):
                    ts, te = (2 if left == "ThreePrime" else 0), len(axis) - (3 if right == "FivePrime" else 0)
                    bs, be = (2 if left == "FivePrime" else 0), len(axis) - (3 if right == "ThreePrime" else 0)
                    value = fragment()
                    value["ends"] = {
                        "left": {"kind": left, "overhang": axis[:2] if left == "FivePrime" else str(Seq(axis[:2]).reverse_complement()), "enzymes": []},
                        "right": {"kind": right, "overhang": str(Seq(axis[-3:]).reverse_complement()) if right == "FivePrime" else axis[-3:], "enzymes": []},
                        "opposite_sequence": str(Seq(axis[bs:be]).reverse_complement()),
                    }
                    value["sequence_qualities"] = value["opposite_sequence_qualities"] = None
                    data = fixture(payload=value, record_edit=lambda r: setattr(r, "seq", Seq(axis[ts:te])))
                    report, _, _ = checker.inspect(data)
                    self.assertEqual(report["geometry"], {"extent": len(axis), "top": {"start": ts, "end": te},
                                                          "bottom": {"start": bs, "end": be}})


if __name__ == "__main__":
    unittest.main()
