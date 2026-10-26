"""Focused standard-library tests for the independent carrier shape checker.

Run: python3 -m unittest discover -s validation -p test_annotation_carrier.py

Synthetic envelopes below have valid framing and record hashes, not recomputed
biological/context digests. They test the checker's expressly limited scope.
"""
import copy
import hashlib
import json
import unittest
import xml.etree.ElementTree as ET

import check_annotation_carrier as checker


UUID = "07e4fc98-b35b-4d51-b9d4-22faf05f5503"


def fixture():
    feature = {
        "id": 1,
        "label": "annotation",
        "kind": "misc_feature",
        "parts": [{"start": 0, "end": 4}],
        "location": {
            "Range": {"start": 0, "end": 4, "before": False, "after": False}
        },
        "reverse": False,
        "qualifiers": [["label", "annotation"]],
        "edited": False,
    }
    reference = dict.fromkeys(checker.REFERENCE_KEYS, None)
    reference.update(title="Original citation", extra=[["DOI", "example"]])
    return {
        "version": 1,
        "instance_id": UUID,
        "tracking_limited": False,
        "feature_epochs": [[1, UUID]],
        "origins": [{
            "id": "1" * 64,
            "source_id": UUID,
            "occurrence_id": "2" * 64,
            "top": {"start": 0, "end": 4},
            "bottom": {"start": 0, "end": 4},
            "feature": copy.deepcopy(feature),
            "references": [reference],
            "extent": 4,
            "molecule": "Dna",
            "circular": False,
            "source_digest": "3" * 64,
            "sequence_digest": "4" * 64,
        }],
        "pieces": [{
            "mapping_digest": "5" * 64,
            "origin_id": "1" * 64,
            "copy_id": UUID,
            "feature_id": 1,
            "axis": "Top",
            "runs": [{"leaf": 0, "offset": 0, "length": 4,
                      "target_start": 0, "reverse": False, "strands": 3}],
            "baseline": feature,
            "references_digest": "6" * 64,
            "sequence_digest": "7" * 64,
        }],
    }


class MemoryFile:
    def __init__(self, data):
        self.data = data

    def read_bytes(self):
        return self.data

    def __str__(self):
        return "<synthetic in-memory carrier>"


def native(registry):
    def packet(kind, body):
        return bytes([kind]) + len(body).to_bytes(4, "big") + body

    features = b'<Features><Feature name="annotation" type="misc_feature"/></Features>'
    records = b"\x0a" + len(features).to_bytes(8, "big") + features
    value = {
        "format": "native",
        "records_sha256": hashlib.sha256(records).hexdigest(),
        "context_sha256": "0" * 64,
        "registry": registry,
    }
    text = json.dumps(value, separators=(",", ":"), ensure_ascii=True)
    text = "".join(f"\\u{ord(c):04x}" if c.isspace() else c for c in text)
    raw = text.encode()
    header = f"bytes={len(raw)};sha256={hashlib.sha256(raw).hexdigest()};json="
    notes = ET.Element("Notes")
    ET.SubElement(notes, "Comments").text = (
        "[PopGene lineage v1 begin]\n" + header + text + "\n[PopGene lineage end]"
    )
    return MemoryFile(
        packet(9, b"SnapGene\x00\x01\x00\x0e\x00\x01")
        + packet(0, b"\x01ACGT") + packet(10, features)
        + packet(6, ET.tostring(notes, encoding="utf-8"))
    )


class CarrierShapeTests(unittest.TestCase):
    def test_valid_nested_shape_reports_only_the_promised_scope(self):
        report = checker.read(native(fixture()))
        self.assertTrue(report["recursive_typed_registry_shape_verified"])
        self.assertFalse(report["ancestor_document_or_sequence_fields_present"])
        self.assertFalse(report["semantic_context_digest_verified"])
        self.assertFalse(report["restoration_semantics_verified"])

    def test_forged_nested_objects_fail_even_with_matching_checksums_and_records(self):
        mutations = {
            "origin object": lambda r: r["origins"][0].update(Document={"sequence": "ACGT"}),
            "location variant": lambda r: r["origins"][0]["feature"].update(
                location={"Document": {"sequence": "ACGT"}}),
            "nested location field": lambda r: r["pieces"][0]["baseline"].update(
                location={"Complement": {"Join": [{"Range": {
                    "start": 0, "end": 4, "before": False, "after": False,
                    "sequence": "ACGT"}}]}}),
            "reference field": lambda r: r["origins"][0]["references"][0].update(
                Document={"sequence": "ACGT"}),
            "reference scalar": lambda r: r["origins"][0]["references"][0].update(
                title={"sequence": "ACGT"}),
            "reference extra pair": lambda r: r["origins"][0]["references"][0].update(
                extra=[["parent", {"sequence": "ACGT"}]]),
            "segment field": lambda r: r["pieces"][0]["baseline"]["parts"][0].update(
                Document={"sequence": "ACGT"}),
            "geometry scalar": lambda r: r["origins"][0]["top"].update(
                end={"sequence": "ACGT"}),
            "qualifier object": lambda r: r["pieces"][0]["baseline"]["qualifiers"].append(
                ["parent", {"sequence": "ACGT"}]),
            "run field": lambda r: r["pieces"][0]["runs"][0].update(sequence="ACGT"),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                registry = fixture()
                mutate(registry)
                with self.assertRaises(AssertionError):
                    checker.read(native(registry))

    def test_scalar_types_and_native_record_indices_are_checked(self):
        mutations = {
            "bool is not version": lambda r: r.update(version=True),
            "bool is not coordinate": lambda r: r["origins"][0]["top"].update(start=False),
            "integer is not strand bool": lambda r: r["pieces"][0]["runs"][0].update(reverse=0),
            "leaf overflow": lambda r: r["pieces"][0]["runs"][0].update(leaf=2**32),
            "negative offset": lambda r: r["pieces"][0]["runs"][0].update(offset=-1),
            "missing optional field": lambda r: r["origins"][0]["references"][0].pop("authors"),
            "noncanonical null location": lambda r: r["pieces"][0]["baseline"].update(location=None),
            "zero native index": lambda r: (r["pieces"][0].update(feature_id=0),
                                            r["pieces"][0]["baseline"].update(id=0)),
            "past-end native index": lambda r: (r["pieces"][0].update(feature_id=2),
                                                r["pieces"][0]["baseline"].update(id=2)),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                registry = fixture()
                mutate(registry)
                with self.assertRaises(AssertionError):
                    checker.read(native(registry))

    def test_all_documented_location_variants_have_typed_nested_payloads(self):
        leaf = {"Range": {"start": 0, "end": 4, "before": True, "after": False}}
        variants = [leaf, {"Between": {"left": 1, "right": 2}},
                    {"Complement": leaf},
                    *[{tag: [leaf]} for tag in ("Join", "Order", "Bond", "OneOf")],
                    {"Remote": {"accession": "AB123.1", "location": leaf}},
                    {"Remote": {"accession": "AB123.1", "location": None}},
                    {"Gap": {"length": 4, "unknown100": False}},
                    {"Gap": {"length": None, "unknown100": True}}]
        for value in variants:
            with self.subTest(value=value):
                checker.location_shape(value)
        for tag in ("Complement", "Join", "Order", "Bond", "OneOf", "Remote"):
            bad = {"Document": {"sequence": "ACGT"}}
            wrapped = {tag: bad if tag == "Complement" else
                       {"accession": "AB123.1", "location": bad} if tag == "Remote" else [bad]}
            with self.subTest(invalid=tag), self.assertRaises(AssertionError):
                checker.location_shape(wrapped)


if __name__ == "__main__":
    unittest.main()
