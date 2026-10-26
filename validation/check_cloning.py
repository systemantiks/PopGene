#!/usr/bin/env python3
"""Compare actual Rust cloning products with independent Bio.Seq/Bio.Restriction
sequence constructions. This validates sequence/coordinate semantics, not wet-lab
yield, mismatch repair, reaction kinetics, or proprietary application parity.
"""
import json
import subprocess
from pathlib import Path
import Bio
from rust_probe import build_example
from Bio.Seq import Seq
from Bio.Restriction import BsaI, EcoRI

ROOT = Path(__file__).resolve().parents[1]

def main():
    assert Bio.__version__ == "1.88"
    probe = build_example('popgene-core', 'cloning_probe')
    cases, expected = [], []
    sequence = Seq((ROOT / "validation/fixtures/primer-design-template.txt").read_text().strip())
    for offset in [0, 225, 150]:
        rotated = sequence[offset:] + sequence[:offset]
        start = (139 - offset) % len(sequence)
        rstart = (394 - offset) % len(sequence)
        f, r = "GGGCGG" + str(sequence[139:159]), "TTAACC" + str(sequence[394:414].reverse_complement())
        cases.append(dict(method="pcr", template=dict(sequence=str(rotated), circular=True), forward=dict(start=start,length=20,sequence=f), reverse=dict(start=rstart,length=20,sequence=r)))
        expected.append(("GGGCGG" + str(sequence[139:414]) + "GGTTAA", False, {"forward":f,"reverse":r}))
    substrate = Seq("AACGGAATTCTGGA")
    assert EcoRI.search(substrate) == [6]
    assert tuple(map(str, EcoRI.catalyse(substrate))) == (str(substrate[:5]), str(substrate[5:]))
    for reverse in [False, True]:
        cases.append(dict(method="digest_ligation",template=dict(sequence=str(substrate)),enzyme="EcoRI",fragments=[1,0] if reverse else [0,1],reverse=reverse))
        expected.append((str(substrate.reverse_complement() if reverse else substrate),False,{}))
    ab, ba = "ACGTCAGTGCATCGATGCTA", "TGCACTAGGCTAACGTCGAT"
    a, b = ba + "CCGAGTAC" + ab, ab + "GGTCATCA" + ba
    for circular in [False,True]:
        cases.append(dict(method="overlap",inputs=[dict(sequence=a),dict(sequence=b,annotations=[dict(label="closing",start=len(b)-25,end=len(b))])],minimum=15,maximum=25,circular=circular))
        expected.append((a + b[len(ab):len(b)-len(ba) if circular else len(b)],circular,{"closing":b[-25:]}))
    cassette=lambda left,body,right:"AAAAGGTCTCA"+left+body+right+"TGAGACCAAAA"
    a,b=cassette("AGTC","CCGATATACCGG","CGTA"),cassette("CGTA","TTAAGGCCAT","AGTC")
    fragments=[BsaI.catalyse(Seq(s)) for s in [a,b]]
    assembled=fragments[0][1]+fragments[1][1]
    assert not BsaI.search(assembled,linear=False)
    cases.append(dict(method="golden_gate",inputs=[dict(sequence=a),dict(sequence=b)],fragments=[1,1],enzyme="BsaI"))
    expected.append((str(assembled),True,{}))
    rows=subprocess.run([str(probe)],input="".join(json.dumps(v)+"\n" for v in cases),capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(rows)==len(cases)
    for i,(row,(sequence,circular,annotations)) in enumerate(zip(rows,expected)):
        got=json.loads(row)
        assert "error" not in got,(i,got)
        assert got["sequence"]==sequence,(i,"sequence",got["sequence"],sequence)
        assert got["circular"]==circular,(i,"topology")
        if got["opposite"] is not None:
            assert got["opposite"]==str(Seq(sequence).reverse_complement()),(i,"opposite")
        for label,expected_sequence in annotations.items():
            feature=next(f for f in got["annotations"] if f["label"]==label)
            extracted=Seq("".join(sequence[p["start"]:p["end"]] for p in feature["parts"]))
            if feature["reverse"]:extracted=extracted.reverse_complement()
            assert str(extracted)==expected_sequence==feature["sequence"],(i,label,feature)
    result=dict(status="pass",biopython=Bio.__version__,production_cases=len(cases),scope=__doc__.strip(),primary_references=["https://www.neb.com/tools-and-resources/usage-guidelines/cloning-guide","https://www.neb.com/en-ca/tools-and-resources/selection-charts/compatible-cohesive-ends-and-generation-of-new-restriction-sites","https://intl.neb.com/en-us/-/media/nebus/files/manuals/manuale1602.pdf","https://biopython.org/docs/1.88/api/Bio.Restriction.Restriction.html"])
    (ROOT/"validation/results").mkdir(exist_ok=True)
    (ROOT/"validation/results/cloning-oracle.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
