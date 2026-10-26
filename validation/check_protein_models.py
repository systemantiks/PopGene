#!/usr/bin/env python3
"""Compare actual Rust protein models with independent Biopython 1.88 calculations.

Run: .venv/bin/python validation/check_protein_models.py
Published EMBOSS and Grimsley/IPC pKa parameters are inserted into Biopython's
independent charge implementation; neither program is invoked as an executable.
Grimsley means come from doi:10.1002/pro.19, with Arg 12.04 supplied by Table 4
of doi:10.1186/s13062-016-0159-9, not measured in the Grimsley summary. The values
are not fitted to a viewer result. Masses use Biopython and Unimod facts. A pI
requires a strict charge sign change within pH 0–14; selected intervals have new
free termini before declared modifications are applied. Unknown residues/stops
and unsupported phosphorylation charge/optics remain unavailable. Secondary
checks cover simple propensity cases, not empirical prediction accuracy.
"""
import json
from pathlib import Path
import random
import subprocess

import Bio
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from rust_probe import build_example

ROOT = Path(__file__).resolve().parents[1]
STANDARD = "ACDEFGHIKLMNPQRSTVWY"
MODELS = ("bjellqvist", "emboss", "grimsley2009")
# Factual parameters transcribed from primary source tables, not production code.
GRIMSLEY_PKA = {"Nterm": 7.7, "Cterm": 3.3, "C": 6.8, "D": 3.5, "E": 4.2,
                "H": 6.6, "K": 10.5, "R": 12.04, "Y": 10.3}
DELTAS = {"n": (42.0367, 42.010565), "k": (42.0367, 42.010565),
          "c": (-0.9848, -0.984016), "m": (15.9994, 15.994915),
          "p": (79.9799, 79.966331), "d": (-2.0159, -2.015650)}


def charge_oracle(selected, model, modifications):
    result = IsoelectricPoint(selected)
    if model == "emboss":
        result.pos_pKs = dict(Nterm=7.5, K=10.8, R=12.5, H=6.5)
        result.neg_pKs = dict(Cterm=3.6, D=3.9, E=4.1, C=8.5, Y=10.1)
    elif model == "grimsley2009":
        result.pos_pKs = {key: GRIMSLEY_PKA[key] for key in ("Nterm", "K", "R", "H")}
        result.neg_pKs = {key: GRIMSLEY_PKA[key] for key in ("Cterm", "D", "E", "C", "Y")}
    for mod in modifications:
        if mod == "n":
            result.charged_aas_content["Nterm"] = 0
        elif mod == "c":
            result.charged_aas_content["Cterm"] = 0
        elif mod.startswith("k:"):
            result.charged_aas_content["K"] -= 1
        elif mod.startswith("d:"):
            result.charged_aas_content["C"] -= 2
    return result


def high_precision_root(charge):
    """Independently bracket the monotone HH charge curve; no boundary clamping."""
    low, high = 0., 14.
    if not (charge(low) > 0 and charge(high) < 0):
        return None
    for _ in range(80):
        mid = (low + high) / 2
        if charge(mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def main():
    assert Bio.__version__ == "1.88", "Use the pinned Biopython 1.88 environment"
    assert __debug__, "Do not run this assertion-based checker with python -O"
    probe = build_example("popgene-core", "protein_probe")
    canonical = "".join((ROOT / "validation/fixtures/regression/protein-canonical.faa")
                        .read_text().splitlines()[1:])
    rng = random.Random(261025)
    sequences = [canonical, "G", "GG", "MCCY", "D" * 100, "K" * 100,
                 "ACDEFGHIKLMNPQRSTVWYKK"]
    sequences += [aa for aa in "CDEHKRY"]
    sequences += ["".join(rng.choices(STANDARD, k=n)) for n in (3, 20, 50, 250)]
    cases = [(seq, 0, len(seq), model, ph, "", 0)
             for seq in sequences for model in MODELS for ph in (0, 2, 7, 9, 14)]
    for model in MODELS:
        cases.extend([
            ("KGGD", 1, 3, model, 7, "", 0),
            ("MCCY", 0, 4, model, 7, "n,c,m:0,d:1:2", 0),
            ("AKAAK", 0, 5, model, 7, "k:1,k:4", 0),
            ("AKCCMYRK", 1, 7, model, 7, "n,c,k:1,d:2:3,m:4", 0),
            (canonical, 0, len(canonical), model, 7, "n,c,k:8,d:1:21", 0),
            ("SY", 0, 2, model, 7, "p:0", 0),
            ("GG", 0, 2, model, 7, "n,c", 0),
            ("DG", 0, 2, model, 7, "n,c", 0),
            ("KG", 0, 2, model, 7, "n,c", 0),
        ])
        cases.extend((seq, 0, len(seq), model, 7, "", 0)
                     for seq in ("AX", "AB", "AZ", "AJ", "AU", "AO", "A*", "A*A"))
    cases.extend((seq, 0, len(seq), "bjellqvist", 7, "", 1)
                 for seq in ("AAAAAAAAAA", "VVVVVVVVVV", "GGGGGGGGGG", "AX*"))
    output = subprocess.run(
        [str(probe), "--detailed"],
        input="".join("\t".join(map(str, case)) + "\n" for case in cases),
        text=True, capture_output=True, check=True,
    ).stdout.splitlines()
    assert len(output) == len(cases)
    canonical_results = {}
    for case, row in zip(cases, output):
        seq, start, end, model, ph, mods, secondary = case
        selected, values = seq[start:end], row.split("\t")
        if not set(selected) <= set(STANDARD):
            assert values == ["-"] * 7, case
            continue
        modifications = mods.split(",") if mods else []
        avg = molecular_weight(selected, "protein") + sum(DELTAS[m[0]][0] for m in modifications)
        mono = molecular_weight(selected, "protein", monoisotopic=True) + sum(
            DELTAS[m[0]][1] for m in modifications)
        assert abs(float(values[0]) - avg) < 1e-6, (case, "average", values[0], avg)
        assert abs(float(values[1]) - mono) < 1e-6, (case, "mono", values[1], mono)
        if any(mod.startswith("p:") for mod in modifications):
            assert values[2:6] == ["-"] * 4, case
            continue
        oracle = charge_oracle(selected, model, modifications)
        charge = oracle.charge_at_pH(ph)
        assert abs(float(values[3]) - charge) < 1e-8, (case, "charge", values[3], charge)
        pi = high_precision_root(oracle.charge_at_pH)
        if pi is None:
            assert values[2] == "-", case
        else:
            assert abs(float(values[2]) - pi) < 1e-9, (case, "pI", values[2], pi)
            assert abs(float(values[2]) - oracle.pi(min_=0, max_=14)) < 1e-4, case
        extinction = (selected.count("W") * 5500 + selected.count("Y") * 1490
                      + sum(mod.startswith("d:") for mod in modifications) * 125)
        assert float(values[4]) == extinction, case
        assert abs(float(values[5]) - extinction / avg) < 1e-8, case
        if secondary:
            assert values[6] == {"A": "H", "V": "E", "G": "C"}[selected[0]] * len(selected)
        if seq == canonical and not mods and ph == 7:
            canonical_results[model] = dict(
                average_da=float(values[0]), pI=float(values[2]), charge_pH7=float(values[3]),
                extinction_280=float(values[4]), absorbance_1mg_ml=float(values[5]),
                independent_pI=pi, independent_charge_pH7=charge)
    assert f"{canonical_results['grimsley2009']['pI']:.2f}" == "7.60"
    assert f"{canonical_results['grimsley2009']['charge_pH7']:.2f}" == "1.18"
    result = dict(
        status="pass", biopython=Bio.__version__, production_cases=len(cases),
        cases_per_model={model: sum(case[3] == model for case in cases) for model in MODELS},
        scope=__doc__.strip(), grimsley_2009_ipc_pka=GRIMSLEY_PKA,
        canonical_results=canonical_results,
        primary_references=[
            "https://web.expasy.org/compute_pi/pi_tool-doc.html",
            "https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat",
            "https://doi.org/10.1002/pro.19", "https://doi.org/10.1186/s13062-016-0159-9",
        ] + [f"https://www.unimod.org/modifications_view.php?editid1={i}"
             for i in (1, 2, 21, 35, 2020)] + [
            "https://web.expasy.org/protscale/pscale/alpha-helixFasman.html",
            "https://web.expasy.org/protscale/pscale/beta-sheetFasman.html"],
    )
    path = ROOT / "validation/results/protein-models-oracle.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
