#!/usr/bin/env python3
"""Production RustBio alignments against independent Biopython1.88 scores.
Global/local/query-in-reference, both query strands, exact circular rotations;
all returned residue coordinates and reconstructed affine path scores checked.
"""
from pathlib import Path
import Bio,json,random,subprocess,tempfile
from rust_probe import compile_probe
from Bio.Align import PairwiseAligner,substitution_matrices
from Bio.Seq import Seq
ROOT=Path(__file__).resolve().parents[1]
RUST=r'''
use popgene_core::{Document,Molecule,alignment::{align_pair,AlignmentMode,AlignmentOptions}};
use std::io::{self,BufRead};
fn main(){for line in io::stdin().lock().lines(){let line=line.unwrap();let p:Vec<_>=line.split_whitespace().collect();let molecule=if p[3]=="1"{Molecule::Protein}else{Molecule::Dna};let mut reference=Document::new_typed("r",p[4],molecule).unwrap();reference.circular=p[2]=="1";let mut query=Document::new_typed("q",p[5],molecule).unwrap();query.sequence_qualities=Some((0..query.sequence.len()).map(|i|Some((i%94) as u8)).collect());let options=AlignmentOptions{mode:match p[0]{"0"=>AlignmentMode::Global,"1"=>AlignmentMode::Local,_=>AlignmentMode::QueryInReference},both_strands:p[1]=="1",..Default::default()};let a=align_pair(&reference,&query,options).unwrap();println!("{}\t{}\t{}\t{}\t{:?}\t{:?}\t{:?}",a.score,a.query.reverse,a.reference.sequence,a.query.sequence,a.reference.positions,a.query.positions,a.query.qualities);}}
'''
def main():
    rng=random.Random(261026);cases=[]
    for i in range(360):
        protein=i%7==0;alphabet='ACDEFGHIKLMNPQRSTVWY' if protein else 'ACGT'
        r=''.join(rng.choices(alphabet,k=rng.randrange(2,28)));q=''.join(rng.choices(alphabet,k=rng.randrange(1,24)))
        cases.append((i%3,i%2,int(i%5==0 and not protein),int(protein),r,q))
    cases += [(2,1,1,0,'AAAACCCC','CCAAAA'),(2,1,1,0,'ATAT','ATATATAT'),(0,0,0,0,'AAACCC','AAA'),(1,0,0,1,'HEAGAWGHEE','PAWHEAE')]
    with tempfile.TemporaryDirectory(prefix='popgene-alignment-oracle-') as directory:
        directory=Path(directory);source=directory/'probe.rs';binary=directory/'probe';source.write_text(RUST)
        compile_probe(source,binary,['popgene-core'],optimize=True)
        output=subprocess.run([str(binary)],input='\n'.join(' '.join(map(str,c)) for c in cases)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(cases)
    matrix=substitution_matrices.load('BLOSUM62')
    for case,line in zip(cases,output):
        mode,both,circular,protein,reference,query=case
        score,reverse,rs,qs,rp,qp,quality=line.split('\t');score=int(score);reverse=reverse=='true'
        decode=lambda text:json.loads(text.replace('Some(','').replace(')','').replace('None','null'))
        rp,qp,quality=map(decode,(rp,qp,quality))
        aligner=PairwiseAligner();aligner.mode='local' if mode==1 else 'global';aligner.open_gap_score=-5;aligner.extend_gap_score=-1
        if protein:aligner.substitution_matrix=matrix
        else:aligner.match_score=2;aligner.mismatch_score=-3
        if mode==2:
            aligner.open_left_deletion_score=0;aligner.extend_left_deletion_score=0
            aligner.open_right_deletion_score=0;aligner.extend_right_deletion_score=0
        references=[reference[i:]+reference[:i] for i in range(len(reference))] if circular and mode!=0 else [reference]
        queries=[query,str(Seq(query).reverse_complement())] if both and not protein else [query]
        expected=max(aligner.score(r,q) for r in references for q in queries)
        assert score==expected,(case,score,expected,rs,qs)
        assert len(rs)==len(qs)==len(rp)==len(qp)==len(quality)
        consumed_r=[];consumed_q=[];path_score=0;gap=None
        for a,b,ri,qi,qual in zip(rs,qs,rp,qp,quality):
            assert (a=='-')==(ri is None);assert (b=='-')==(qi is None)
            if ri is not None:assert a==reference[ri];consumed_r.append(ri)
            if qi is not None:
                expected_base=str(Seq(query[qi]).complement()) if reverse else query[qi]
                assert b==expected_base;assert qual==qi%94;consumed_q.append(qi)
            else:assert qual is None
            next_gap='r' if a=='-' else 'q' if b=='-' else None
            if next_gap:path_score+=-1 if gap==next_gap else -5
            else:path_score+=int(matrix[a,b]) if protein else 2 if a==b else -3
            gap=next_gap
        assert path_score==score,(case,path_score,score)
        assert len(consumed_r)==len(set(consumed_r)),('more than one circle lap',case)
        assert len(consumed_q)==len(set(consumed_q))
        if mode!=1:assert sorted(consumed_q)==list(range(len(query)))
        if mode==0:assert consumed_r==list(range(len(reference)))
    result={'status':'pass','biopython':Bio.__version__,'production_cases':len(cases),'scope':__doc__.strip(),'reference':'https://biopython.org/docs/1.88/Tutorial/chapter_pairwise.html','protein_matrix':'BLOSUM62','gap_formula':'open + (length - 1) * extend'}
    destination=ROOT/'validation/results/alignment-oracle.json';destination.parent.mkdir(exist_ok=True);destination.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
