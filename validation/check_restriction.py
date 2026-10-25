#!/usr/bin/env python3
"""Compare the actual Rust restriction engine with Biopython's independent engine.

Run with .venv/bin/python validation/check_restriction.py. No network is used.
The standalone Rust probe imports the production module, so this does not need
the desktop dependency graph or a separately maintained implementation.
"""

import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

import Bio
from Bio import Restriction, SeqIO
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "crates/popgene-core/src/restriction.rs"
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
    "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
    "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def find_rustc():
    configured = os.environ.get("RUSTC")
    if configured:
        return configured
    local = list((ROOT / ".tools/rustup/toolchains").glob("stable-*/bin/rustc"))
    return str(local[0]) if local else shutil.which("rustc") or "rustc"


def main():
    rng = random.Random(823)
    source = f'#[path = {json.dumps(str(MODULE))}]\nmod restriction;\n' + r'''
use std::io::{self, BufRead};
fn main() {
    for enzyme in restriction::ENZYMES {
        println!("E\t{}\t{}\t{}\t{}", enzyme.name, enzyme.recognition,
            enzyme.top_cut, enzyme.bottom_cut);
    }
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let parts: Vec<_> = line.splitn(3, '\t').collect();
        let start = std::time::Instant::now();
        let result = restriction::find_sites(parts[2], parts[1] == "circular");
        println!("T\t{}\t{}\t{}", parts[0], start.elapsed().as_micros(), result.truncated);
        for (index, count) in result.counts.iter().enumerate() {
            println!("C\t{}\t{}\t{}", parts[0], restriction::ENZYMES[index].name, count);
        }
        for site in result.sites {
            let top = site.top_cut.map_or("-".into(), |n| n.to_string());
            let bottom = site.bottom_cut.map_or("-".into(), |n| n.to_string());
            println!("S\t{}\t{}\t{}\t{}\t{}\t{}", parts[0],
                restriction::ENZYMES[site.enzyme_index].name,
                site.start, site.reverse, top, bottom);
        }
    }
}
'''
    with tempfile.TemporaryDirectory(prefix="popgene-restriction-") as temporary:
        temporary = Path(temporary)
        (temporary / "probe.rs").write_text(source)
        probe = temporary / "probe"
        subprocess.run([find_rustc(), "--edition=2024", "-O", "-o", str(probe),
                        str(temporary / "probe.rs")], check=True)
        catalogue_output = subprocess.run([str(probe)], input="", text=True,
                                          capture_output=True, check=True).stdout
        catalogue = []
        for line in catalogue_output.splitlines():
            _, name, motif, top, bottom = line.split("\t")
            catalogue.append((name, motif, int(top), int(bottom)))

        missing = []
        for name, motif, top, bottom in catalogue:
            enzyme = getattr(Restriction, name, None)
            if enzyme is None:
                missing.append(name)
                continue
            assert (motif, top, bottom) == (
                enzyme.site, enzyme.fst5, len(enzyme.site) + enzyme.fst3
            ), (name, "catalogue differs from Biopython")
        assert missing == ["EcoP15I"], missing

        def random_dna(length):
            return "".join(rng.choices("ACGT", k=length))

        cases = []
        fixture = str(SeqIO.read(ROOT / "validation/fixtures/annotated.gb", "genbank").seq)
        for filename in ("annotated.dna",):
            roundtrip = SeqIO.read(ROOT / "validation/fixtures" / filename, "snapgene")
            assert str(roundtrip.seq).upper() == fixture.upper(), filename
        cases.append(("annotated-linear", False, fixture))
        cases.append(("annotated-circular", True, fixture))
        random_sequence = random_dna(40_000)
        cases.append(("random-linear", False, random_sequence))
        cases.append(("random-circular", True, random_sequence))
        cases.append(("lowercase", False, random_sequence.lower()))
        cases.append(("ambiguous-targets", False, "N".join(
            random_sequence[i:i + 19] for i in range(0, 4000, 19))))
        for name, motif, _, _ in catalogue:
            for variant in range(4):
                concrete = "".join(rng.choice(IUPAC[base]) for base in motif)
                orientation = concrete if variant % 2 == 0 else str(Seq(concrete).reverse_complement())
                sequence = random_dna(80) + orientation + random_dna(80)
                cases.append((f"{name}-{variant}-linear", False, sequence))
                # Place the chosen recognition site across the circular origin.
                rotation = 80 + len(orientation) // 2
                sequence = sequence[rotation:] + sequence[:rotation]
                cases.append((f"{name}-{variant}-circular", True, sequence))

        rows = "".join(f"{name}\t{'circular' if circular else 'linear'}\t{seq}\n"
                       for name, circular, seq in cases)
        output = subprocess.run([str(probe)], input=rows, text=True,
                                capture_output=True, check=True).stdout
        counts, sites, timings = {}, {}, {}
        for line in output.splitlines():
            fields = line.split("\t")
            if fields[0] == "C":
                counts[(fields[1], fields[2])] = int(fields[3])
            elif fields[0] == "S":
                case, name, start, reverse, top, bottom = fields[1:]
                sites.setdefault((case, name), []).append((
                    int(start), reverse == "true",
                    None if top == "-" else int(top),
                    None if bottom == "-" else int(bottom)))
            elif fields[0] == "T":
                timings[fields[1]] = int(fields[2])
                assert fields[3] == "false", (fields[1], "unexpected truncation")

        comparisons = 0
        ambiguous_oracle_hits_excluded = 0
        for case, circular, sequence in cases:
            n = len(sequence)
            for name, _, _, _ in catalogue:
                if name in missing:
                    continue
                enzyme = getattr(Restriction, name)
                actual = sites.get((case, name), [])
                assert counts[(case, name)] == len(actual)
                # Biopython omits cuts at terminal boundaries 0 and n as well
                # as cuts beyond the molecule. PopGene retains recognition
                # sites there; compare only Biopython's supported interior.
                actual = sorted((top, bottom) for _, _, top, bottom in actual
                                if top is not None and bottom is not None and
                                (circular or (0 < top < n and 0 < bottom < n)))
                expected = []
                for after_cut in enzyme.search(Seq(sequence), linear=not circular):
                    top = after_cut - 1
                    bottom = top - enzyme.ovhg
                    if circular:
                        top, bottom = top % n, bottom % n
                    expected.append((top, bottom))
                if any(base.upper() not in "ACGT" for base in sequence):
                    # Biopython's regex uses '.' for motif N and can therefore
                    # match an unknown target N. Its independent regex supplies
                    # the recognition spans; keep only definite A/C/G/T spans.
                    scan = sequence.upper()
                    if circular:
                        scan += scan[:len(enzyme.site) - 1]
                    definite_top_cuts = set()
                    for match in enzyme.compsite.finditer(scan):
                        forward = match.groupdict().get(name)
                        recognized = forward or match.groupdict().get(name + "_as")
                        if recognized and all(base in "ACGT" for base in recognized):
                            top = match.start() + (enzyme.fst5 if forward else -enzyme.fst3)
                            definite_top_cuts.add(top % n if circular else top)
                    filtered = [pair for pair in expected if pair[0] in definite_top_cuts]
                    ambiguous_oracle_hits_excluded += len(expected) - len(filtered)
                    expected = filtered
                assert actual == sorted(expected), (case, name, actual, expected)
                comparisons += 1

        ecori = sites[("annotated-circular", "EcoRI")]
        assert ecori == [(2344, False, 2345, 2349)], ecori
        nmeaiii = sites[("annotated-circular", "NmeAIII")]
        assert nmeaiii == [(1602, True, 1583, 1581)], nmeaiii
        # EcoP15I is absent from this Biopython release; its catalogue and
        # orientation/boundary values have separate primary-source unit tests.
        large = str(SeqIO.read(ROOT / "validation/fixtures/large-5000000.fasta", "fasta").seq)
        benchmark_output = subprocess.run([str(probe)],
            input=f"large\tlinear\t{large}\n", text=True, capture_output=True, check=True).stdout
        benchmark = next(line.split("\t") for line in benchmark_output.splitlines()
                         if line.startswith("T\t"))
        result = {
            "biopython_version": Bio.__version__,
            "catalogue_size": len(catalogue),
            "catalogue_entries_compared": len(catalogue) - len(missing),
            "not_in_biopython": missing,
            "synthetic_cases": len(cases),
            "enzyme_case_comparisons": comparisons,
            "ambiguous_biopython_hits_excluded": ambiguous_oracle_hits_excluded,
            "fixture_ecori": {"recognition_start_zero_based": 2344,
                              "top_cut_boundary_zero_based": 2345,
                              "biopython_next_base_one_based": 2346},
            "fixture_nmeaiii": {"recognition_start_zero_based": 1602,
                                "reverse": True,
                                "top_cut_boundary_zero_based": 1583,
                                "bottom_cut_boundary_zero_based": 1581,
                                "biopython_next_base_one_based": 1584,
                                "nominal_outside_motif_offsets": "21/19",
                                "native_fixture_sequence_equal_ignoring_case": True},
            "large_fixture_bases": len(large),
            "large_fixture_scan_ms": int(benchmark[2]) / 1000,
            "large_fixture_positions_truncated": benchmark[3] == "true",
            "result": "pass",
        }
        result_path = ROOT / "validation/results/restriction-oracle.json"
        result_path.parent.mkdir(parents=True, exist_ok=True)
        result_path.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
