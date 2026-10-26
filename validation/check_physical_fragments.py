"""Independent, read-only checker for the public physical-fragment v1 supplement.

Use the pinned environment (Biopython 1.88 in validation/requirements.lock):
  .venv/bin/python validation/check_physical_fragments.py FILE.gb FILE.embl
  .venv/bin/python validation/check_physical_fragments.py --ecori-matrix DIR

Generate the 18 matrix files in a NEW directory with POPGENE_PHYSICAL_QA_DIR
when running the formats test real_save_reopen_of_cohesive_partial_cds_rejoins_all_three_phases.

Checks framing/checksum, exact nested field shapes, raw FEATURES/FT binding,
ordinary top records, physical strand geometry, metadata bounds and citations.
The optional EcoRI matrix independently checks known sequences, partial CDS
positions/phases and translations using Biopython. This does NOT recompute the
Rust canonical context_sha256, prove ancestry actionability/reconstruction,
validate every location projection/display qualifier, or establish authenticity.
In particular, a valid checksum is not proof that metadata belongs to the top
record; production context validation remains necessary. No Rust is invoked.
"""

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import sys

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, ExactPosition

import check_annotation_carrier as shapes

START = "[PopGene physical fragment v"
END = "[PopGene physical fragment end]"
ASCII_SPACE = " \t\r\n\v\f"
MAX_PAYLOAD = 4 * 1024 * 1024
MAX_FILE = 128 * 1024 * 1024
MAX_BASES = 10_000_000
DNA = frozenset("ACGTRYSWKMBDHVN")
FRAGMENT_KEYS = {
    "ends", "duplex_annotations", "sequence_colors", "sequence_qualities",
    "opposite_sequence_qualities", "methylation", "references",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def object_pairs(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def ordinary_records(data):
    """Extract raw binding bytes and comments without inventing word folding.

    GenBank continuation newlines remain newlines; EMBL CC continuation lines
    are words of the same paragraph. XX separates EMBL comment paragraphs.
    """
    text = data.decode("utf-8")
    if text.startswith("LOCUS"):
        format_name = "genbank"
    elif text.startswith("ID   "):
        format_name = "embl"
    else:
        raise ValueError("expected one GenBank or EMBL record")
    comments, records = [], bytearray()
    in_comment = in_features = False
    for line in text.splitlines(keepends=True):
        if format_name == "genbank":
            if line.startswith("FEATURES"):
                in_features = True
            elif line.startswith(("ORIGIN", "CONTIG")):
                in_features = False
            if in_features:
                records.extend(line.encode())
            if line.startswith("COMMENT     "):
                comments.append(line[12:].strip())
                in_comment = True
            elif in_comment and line.startswith(" " * 12):
                comments[-1] += "\n" + line[12:].strip()
            else:
                in_comment = False
        else:
            if line.startswith("FT"):
                records.extend(line.encode())
            if line.startswith("CC   "):
                if in_comment:
                    comments[-1] += " " + line[5:].strip()
                else:
                    comments.append(line[5:].strip())
                in_comment = True
            else:
                in_comment = False
    record = SeqIO.read(io.StringIO(text), format_name)
    require(record.annotations.get("topology") == "linear", "physical ends require linear topology")
    require(record.annotations.get("molecule_type") == "DNA", "physical ends require double-stranded DNA")
    top = str(record.seq).upper()
    require(top and set(top) <= DNA, "ordinary top is not nonempty IUPAC DNA")
    return format_name, comments, bytes(records), record


def payload(comments, format_name, records):
    starts = [(i, m.start()) for i, c in enumerate(comments) for m in re.finditer(re.escape(START), c)]
    require(len(starts) == 1 and sum(c.count(END) for c in comments) == 1,
            "expected exactly one physical start/end marker pair")
    index, start = starts[0]
    comment = comments[index]
    version = re.match(r"\[PopGene physical fragment v([0-9]+) begin\]", comment[start:])
    require(version is not None and version[1] == "1", "unsupported or malformed physical version")
    begin = start + version.end()
    end = comment.find(END, begin)
    require(end >= begin, "physical markers must be ordered in one comment paragraph")
    framed = comment[begin:end]
    require(not any(c.isspace() and c not in ASCII_SPACE for c in framed), "non-ASCII payload wrapping")
    compact_chars = []
    for character in framed:
        if character not in ASCII_SPACE:
            compact_chars.append(character)
            require(len(compact_chars) <= MAX_PAYLOAD + 128, "physical payload exceeds4MiB")
    compact = "".join(compact_chars)
    match = re.fullmatch(r"bytes=([0-9]{1,10});sha256=([0-9a-f]{64});json=(.*)", compact, re.S)
    require(match is not None, "malformed physical payload header")
    raw = match[3].encode()
    require(len(raw) == int(match[1]) <= MAX_PAYLOAD, "physical payload length/budget mismatch")
    require(sha(raw) == match[2], "physical payload checksum mismatch")
    shapes.depth_check(match[3])
    value = json.loads(raw, object_pairs_hook=object_pairs,
                       parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON number")))
    shapes.exact_keys(value, {"format", "records_sha256", "context_sha256", "fragment"})
    shapes.digest(value["records_sha256"])
    shapes.digest(value["context_sha256"])
    require(value["format"] == format_name, "carrier format differs from ordinary record")
    require(value["records_sha256"] == sha(records), "raw annotation record hash mismatch")
    return value, raw


def segment(value, extent):
    shapes.segment_shape(value)
    require(0 <= value["start"] < value["end"] <= extent, "range outside physical coordinate extent")


def location_bounds(node, extent):
    """Bound local coordinates; remote accessions never become local positions."""
    tag, value = next(iter(node.items()))
    if tag == "Range":
        require(0 <= value["start"] < value["end"] < 2**63, "invalid location range")
        require(extent is None or value["end"] <= extent, "location outside physical extent")
    elif tag == "Between":
        require(value["left"] < 2**63 - 1 and value["right"] < 2**63 - 1, "between endpoint overflow")
        require(value["right"] == value["left"] + 1 or (extent is None and value["right"] == 0),
                "nonadjacent between endpoints on a linear fragment")
        require(extent is None or max(value.values()) < extent, "between site outside physical extent")
    elif tag == "Complement":
        location_bounds(value, extent)
    elif tag in ("Join", "Order", "Bond", "OneOf"):
        require(bool(value), "empty compound location")
        for child in value:
            location_bounds(child, extent)
    elif tag == "Remote":
        require(re.fullmatch(r"[A-Za-z0-9_.]+", value["accession"]), "invalid remote accession")
        if value["location"] is not None:
            location_bounds(value["location"], None)
    elif tag == "Gap":
        require(value["length"] is None or 0 < value["length"] <= MAX_BASES, "gap length outside1..10000000")
        require(not value["unknown100"] or value["length"] is None, "conflicting gap fields")


def citation(value, references):
    require(isinstance(value, str), "citation has no numeric reference slot")
    match = re.fullmatch(r"\[?([0-9]+)\]?", value.strip())
    require(match is not None and 1 <= int(match[1]) <= len(references), "citation outside current references")


def colors(values, extent):
    shapes.array(values, 100_000)
    for value in values:
        shapes.exact_keys(value, {"range", "reverse", "color_rgb"})
        segment(value["range"], extent)
        shapes.boolean(value["reverse"])
        shapes.unsigned(value["color_rgb"], 24)


def qualities(values, length):
    if values is not None:
        shapes.array(values)
        require(len(values) == length, "quality vector length mismatch")
        for value in values:
            require(value is None or (type(value) is int and 0 <= value <= 93),
                    "quality must be null or a Phred score0..93")


def fragment_geometry(fragment, record):
    shapes.exact_keys(fragment, FRAGMENT_KEYS)
    ends = fragment["ends"]
    required = {"left", "right", "opposite_sequence"}
    require(isinstance(ends, dict) and required <= set(ends)
            and set(ends) <= required | {"left_modification", "right_modification"}, "unknown or missing ends field")
    for key in ("left_modification", "right_modification"):
        if key in ends:
            shapes.string(ends[key])  # Canonical None is omitted, not JSON null.
    for side in ("left", "right"):
        end = ends[side]
        shapes.exact_keys(end, {"kind", "overhang", "enzymes"})
        require(end["kind"] in ("Native", "Blunt", "FivePrime", "ThreePrime"), "unknown end polarity")
        shapes.string(end["overhang"])
        require(set(end["overhang"]) <= DNA, "overhang is not uppercase IUPAC DNA")
        require(bool(end["overhang"]) == (end["kind"] in ("FivePrime", "ThreePrime")),
                "overhang length disagrees with end polarity")
        shapes.array(end["enzymes"])
        for enzyme in end["enzymes"]:
            shapes.string(enzyme)
    top = str(record.seq).upper()
    bottom = ends["opposite_sequence"]
    shapes.string(bottom)
    require(set(bottom) <= DNA, "opposite strand is not uppercase IUPAC DNA")
    left, right = ends["left"], ends["right"]
    tp = len(left["overhang"]) if left["kind"] == "FivePrime" else 0
    ts = len(right["overhang"]) if right["kind"] == "ThreePrime" else 0
    bp = len(right["overhang"]) if right["kind"] == "FivePrime" else 0
    bs = len(left["overhang"]) if left["kind"] == "ThreePrime" else 0
    require(tp + ts < len(top) and bp + bs < len(bottom), "duplex has no paired residues")
    require(str(Seq(top[tp:len(top)-ts]).reverse_complement()) == bottom[bp:len(bottom)-bs],
            "paired top/opposite residues disagree")
    require((top[:tp] if tp else bottom[len(bottom)-bs:]) == left["overhang"], "left overhang differs from physical strand")
    require((bottom[:bp] if bp else top[len(top)-ts:]) == right["overhang"], "right overhang differs from physical strand")
    extent = bs + len(top) + bp
    shapes.array(fragment["references"])
    references = fragment["references"]
    require(len(references) == len(record.annotations.get("references", [])), "ordinary/current reference slot count differs")
    for reference in references:
        shapes.reference_shape(reference)
    annotations = fragment["duplex_annotations"]
    shapes.exact_keys(annotations, {"features", "sequence_colors"})
    shapes.array(annotations["features"], 8192)
    for index, feature in enumerate(annotations["features"], len(record.features) + 1):
        shapes.feature_shape(feature)
        require(feature["id"] == index, "supplemental occurrence ID differs from emitted feature order")
        for part in feature["parts"]:
            segment(part, extent)
        if "location" in feature:
            location_bounds(feature["location"], extent)
        for key, value in feature["qualifiers"]:
            if key == "citation":
                citation(value, references)
    for feature in record.features:
        for value in feature.qualifiers.get("citation", []):
            citation(value, references)
    colors(fragment["sequence_colors"], len(top))
    colors(annotations["sequence_colors"], extent)
    qualities(fragment["sequence_qualities"], len(top))
    qualities(fragment["opposite_sequence_qualities"], len(bottom))
    if fragment["methylation"] is not None:
        shapes.exact_keys(fragment["methylation"], {"dam", "dcm", "eco_ki"})
        for value in fragment["methylation"].values():
            shapes.boolean(value)
    return {"extent": extent, "top": {"start": bs, "end": bs + len(top)},
            "bottom": {"start": tp, "end": extent - ts}}


def inspect(data, name="<memory>"):
    require(__debug__, "run without Python -O: shared shape checks use assertions")
    require(Bio.__version__ == "1.88", "use pinned Biopython1.88 from validation/requirements.lock")
    require(len(data) <= MAX_FILE, "file exceeds128MiB")
    format_name, comments, records, record = ordinary_records(data)
    value, raw = payload(comments, format_name, records)
    geometry = fragment_geometry(value["fragment"], record)
    top = str(record.seq).upper()
    opposite = value["fragment"]["ends"]["opposite_sequence"]
    report = {"path": str(name), "sha256": sha(data), "format": format_name,
              "payload_bytes": len(raw), "payload_sha256": sha(raw),
              "records_sha256": sha(records),
              "ordinary_top": top if len(top) <= 4096 else None,
              "ordinary_top_length": len(top), "ordinary_top_sha256": sha(top.encode()),
              "opposite_5to3": opposite if len(opposite) <= 4096 else None,
              "opposite_length": len(opposite), "opposite_sha256": sha(opposite.encode()),
              "geometry": geometry, "supplemental_features": len(value["fragment"]["duplex_annotations"]["features"]),
              "framing_schema_raw_records_and_geometry_verified": True,
              "canonical_context_digest_verified": False, "complete_location_projection_verified": False,
              "ancestry_restoration_verified": False, "authenticity_claimed": False}
    return report, record, value["fragment"]


def ecori_expectation(path, record, fragment):
    match = re.fullmatch(r"ecori-phase-([123])-(fragment-[12]|rejoined)\.(gb|embl)", path.name)
    require(match is not None, "unexpected EcoRI matrix filename")
    phase, part = int(match[1]), match[2]
    source = "ATGGAATTCAAGTAA"
    expected_top = {"fragment-1": "ATGG", "fragment-2": "AATTCAAGTAA", "rejoined": source}[part]
    opposite = {"fragment-1": "AATTCCAT", "fragment-2": "TTACTTG", "rejoined": "TTACTTGAATTCCAT"}[part]
    require(str(record.seq).upper() == expected_top, "unexpected EcoRI product top")
    ends = fragment["ends"]
    require(ends["opposite_sequence"] == opposite, "unexpected EcoRI product opposite")
    for side in ("left", "right"):
        cut = (part == "fragment-1" and side == "right") or (part == "fragment-2" and side == "left")
        require(ends[side] == {"kind": "FivePrime" if cut else "Native", "overhang": "AATT" if cut else "",
                               "enzymes": ["EcoRI"] if cut else []}, "unexpected EcoRI end")
    require(ends.get("left_modification") == ends.get("right_modification") ==
            ("5′ phosphate" if part == "rejoined" else None), "unexpected terminal chemistry")
    annotations = fragment["duplex_annotations"]
    require(annotations["sequence_colors"] == fragment["sequence_colors"] == [], "unexpected fixture colors")
    require(fragment["sequence_qualities"] is None and fragment["opposite_sequence_qualities"] is None,
            "unexpected fixture qualities")
    require(fragment["methylation"] is None and fragment["references"] == [], "unexpected fixture metadata")
    translation = None
    if part == "fragment-1":
        require(not record.features and len(annotations["features"]) == 1, "wrong left fragment annotations")
        feature = annotations["features"][0]
        require(feature["id"] == 1 and feature["kind"] == "CDS" and feature["label"] == "coding α"
                and not feature["reverse"] and feature["edited"], "wrong left CDS identity")
        require(feature["parts"] == [{"start": 0, "end": 8}] and feature["location"] ==
                {"Range": {"start": 0, "end": 8, "before": False, "after": True}}, "wrong left partial CDS")
        require(dict(feature["qualifiers"]) == {"label": "coding α", "partial": None, "codon_start": str(phase)},
                "wrong left phase/qualifiers")
        extraction = expected_top + str(Seq(opposite[:4]).reverse_complement())
        require(extraction == source[:8], "wrong left full-axis coding residues")
    else:
        require(not annotations["features"] and len(record.features) == 1, "wrong ordinary CDS count")
        feature = record.features[0]
        require(feature.type == "CDS" and feature.location.strand == 1
                and int(feature.location.start) == 0 and int(feature.location.end) == len(expected_top), "wrong CDS coordinates")
        require(type(feature.location.start) is (BeforePosition if part == "fragment-2" else ExactPosition)
                and type(feature.location.end) is ExactPosition, "wrong CDS partial endpoints")
        require(int(feature.qualifiers["codon_start"][0]) == ((phase - 1 - 4) % 3 + 1 if part == "fragment-2" else phase),
                "wrong CDS phase")
        extraction = str(feature.extract(record.seq)).upper()
        if part == "rejoined":
            coding = extraction[phase - 1:]
            translation = str(Seq(coding[:len(coding)//3*3]).translate())
            require(translation == {1: "MEFK*", 2: "WNSS", 3: "GIQV"}[phase], "wrong rejoined translation")
    return {"phase": phase, "part": part, "extracted_sequence": extraction, "translation": translation,
            "independent_fixture_expectations_verified": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="*", type=Path)
    parser.add_argument("--ecori-matrix", type=Path, help="also check exactly18 named EcoRI phase/product fixtures")
    parser.add_argument("--output", type=Path, help="write JSON to a new file instead of stdout; existing files are refused")
    args = parser.parse_args()
    paths = list(args.files)
    matrix = set()
    if args.ecori_matrix:
        matrix = {args.ecori_matrix / f"ecori-phase-{phase}-{part}.{ext}"
                  for phase in (1, 2, 3) for part in ("fragment-1", "fragment-2", "rejoined") for ext in ("gb", "embl")}
        require(set(args.ecori_matrix.glob("*.gb")) | set(args.ecori_matrix.glob("*.embl")) == matrix,
                "EcoRI directory must contain exactly the18 named GB/EMBL files")
        paths += sorted(matrix)
    require(bool(paths), "provide FILE arguments or --ecori-matrix DIR")
    reports = []
    for path in paths:
        require(path.stat().st_size <= MAX_FILE, f"{path}: file exceeds128MiB")
        report, record, fragment = inspect(path.read_bytes(), path)
        if path in matrix:
            report["ecori_expectations"] = ecori_expectation(path, record, fragment)
            report["ancestry_shape_and_raw_binding"] = shapes.read(path)
        reports.append(report)
    result = {"status": "pass", "biopython": Bio.__version__, "files": reports,
              "scope": "Framing, exact nested schema, raw annotation binding and physical geometry. Optional known EcoRI fixture expectations. Canonical context/actionable ancestry/authenticity are NOT independently proven.",
              "rust_invoked": False, "input_files_modified": False}
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        with args.output.open("x") as handle:
            handle.write(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, ValueError, OSError, KeyError, TypeError) as error:
        print(f"Physical fragment check failed: {error}", file=sys.stderr)
        raise SystemExit(1) from None
