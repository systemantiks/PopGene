#!/usr/bin/env python3
"""Compare production Rust translations with Biopython and brute-force ORFs.

NCBI's factual code table is independently checked against Biopython 1.88;
Rust's ambiguity/start/stop algorithms are exercised for all 15^3 IUPAC codons
and all 27 codes. ORFs are checked by exhaustive per-start searches, separate
from production's linear interval scanner. No vendor binaries are used.
"""
from pathlib import Path
from rust_probe import compile_probe
import itertools,json,random,subprocess,tempfile
import Bio
from Bio.Data import CodonTable,IUPACData
from Bio.Seq import Seq
ROOT=Path(__file__).resolve().parents[1]
ALPHABET='TCAGYWRKMSBDHVN'

def main():
    reference=json.loads((ROOT/'validation/fixtures/translation/ncbi-genetic-codes.json').read_text())
    ids=[t['id'] for t in reference['tables']]
    assert ids==sorted(CodonTable.unambiguous_dna_by_id)
    concrete=[''.join(c) for c in itertools.product('TCAG',repeat=3)]
    for t in reference['tables']:
        bio=CodonTable.unambiguous_dna_by_id[t['id']]
        assert t['amino_acids']==''.join(bio.forward_table.get(c,'*') for c in concrete)
        assert t['starts']==''.join('M' if c in bio.start_codons else '*' if c in bio.stop_codons else '-' for c in concrete)
    codons=[''.join(c) for c in itertools.product(ALPHABET,repeat=3)]
    source=r'''
use popgene_core::{Document,translation::{GeneticCode,translate_codon,is_start_codon,is_stop_codon,find_orfs,OrfOptions}};
use std::io::{self,BufRead};
fn main(){for line in io::stdin().lock().lines(){let line=line.unwrap();let fields:Vec<_>=line.split_whitespace().collect();
if fields[0]=="T" {let code=GeneticCode::from_id(fields[1].parse().unwrap()).unwrap();let c=fields[2].as_bytes();let c=[c[0],c[1],c[2]];println!("{}{}{}",translate_codon(c,code) as char,if is_start_codon(c,code,true){'M'}else{'-'},if is_stop_codon(c,code){'*'}else{'-'});}
else {let code=GeneticCode::from_id(fields[1].parse().unwrap()).unwrap();let mut doc=Document::new("oracle",fields[5]).unwrap();doc.circular=fields[2]=="1";let options=OrfOptions {genetic_code:code,min_codons:fields[3].parse().unwrap(),alternative_starts:fields[4]=="1",..OrfOptions::default()};let result=find_orfs(&doc,options);println!("{}",result.orfs.iter().map(|o|format!("{}:{}:{}",o.start,u8::from(o.reverse),o.nucleotide_length)).collect::<Vec<_>>().join(","));}}
}
'''
    rng=random.Random(8232026)
    cases=[]
    for i in range(600):
        code=ids[i%len(ids)];circular=i%2;minimum=i%5;alternative=(i//2)%2
        seq=''.join(rng.choices('ACGT',k=rng.randrange(3,100)))
        cases.append((code,circular,minimum,alternative,seq))
    cases += [(1,1,2,0,'AAATAACATG'),(1,0,2,0,'ATGAAATAATTATTTCAT'),(11,0,2,1,'GTGAAATAA'),(1,0,1,0,'ATGTANTAA')]
    rows=[f'T {i} {c}' for i in ids for c in codons]+[f'O {i} {circle} {minimum} {alternative} {seq}' for i,circle,minimum,alternative,seq in cases]
    with tempfile.TemporaryDirectory(prefix='popgene-translation-') as temporary:
        temporary=Path(temporary);rust=temporary/'probe.rs';binary=temporary/'probe';rust.write_text(source)
        compile_probe(rust,binary,['popgene-core'],optimize=True)
        outputs=subprocess.run([str(binary)],input='\n'.join(rows)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(outputs)==len(rows)
    index=0
    for i in ids:
        table=CodonTable.unambiguous_dna_by_id[i]
        for codon in codons:
            possibilities=[''.join(c) for c in itertools.product(*(IUPACData.ambiguous_dna_values[b] for b in codon))]
            residues={table.forward_table.get(c,'*') for c in possibilities}
            expected=(next(iter(residues)) if len(residues)==1 else 'X')+('M' if all(c in table.start_codons for c in possibilities) else '-')+('*' if all(c in table.stop_codons for c in possibilities) else '-')
            assert outputs[index]==expected,(i,codon,outputs[index],expected);index+=1
    for case,output in zip(cases,outputs[index:]):
        code,circular,minimum,alternative,sequence=case;table=CodonTable.unambiguous_dna_by_id[code];n=len(sequence);expected={}
        for reverse,seq in [(False,sequence),(True,str(Seq(sequence).reverse_complement()))]:
            for start in range(n):
                limit=n if circular else n-start
                first=''.join(seq[(start+j)%n] for j in range(3)) if limit>=3 else ''
                starts=table.start_codons if alternative else ['ATG']
                if first not in starts:continue
                for length in range(6,limit+1,3):
                    codon=''.join(seq[(start+length-3+j)%n] for j in range(3))
                    possibilities=[''.join(c) for c in itertools.product(*(IUPACData.ambiguous_dna_values[b] for b in codon))]
                    if all(c in table.stop_codons for c in possibilities):
                        if length//3-1>=minimum:
                            key=(reverse,start%3,start+length)
                            if key not in expected or length>expected[key][2]:expected[key]=(n-1-start if reverse else start,int(reverse),length)
                        break
                    if any(c in table.stop_codons for c in possibilities):break
        actual={tuple(map(int,item.split(':'))) for item in output.split(',') if item}
        assert actual==set(expected.values()),(case,actual,set(expected.values()))
    result={'status':'pass','biopython_version':Bio.__version__,'ncbi_source':reference['source'],'source_sha256':reference['source_sha256'],'active_codes':len(ids),'iupac_codon_cases':index,'independent_orf_cases':len(cases),'orf_method':'brute-force each start to first definite stop, longest per frame/stop interval; one lap on circles'}
    (ROOT/'validation/results').mkdir(exist_ok=True)
    (ROOT/'validation/results/translation-oracle.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
