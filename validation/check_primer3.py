#!/usr/bin/env python3
"""Compare the production C/Rust integration with primer3-py 2.3.1.
The wrappers share Primer3 source lineage, so this validates integration/units,
not an independently derived physical theory. One default Tm also uses Bio's
independent nearest-neighbour implementation with the same published parameters.
"""
import hashlib, json, math, random, subprocess
from rust_probe import build_example
from pathlib import Path
import primer3
from Bio.SeqUtils import MeltingTemp
ROOT=Path(__file__).resolve().parents[1]
def main():
 assert primer3.__version__=='2.3.1'
 probe=build_example('popgene-core','primer_probe')
 rng=random.Random(261025)
 cases=[]
 for model in [0,1]:
  for salt in [0,1,2]:
   for condition in [(50,50,1.5,.6,0,0,37),(250,80,3,.2,2,.5,25)]:
    for sequence in ['CGTTCCAAAGATGTGGGCATGAGCTTAC','GCGCGCTTTAGCGCGC',''.join(rng.choices('ACGT',k=30))]:
     cases.append((str(len(cases)),sequence,*condition,model,salt))
 rows=''.join('T\t'+'\t'.join(map(str,c))+'\n' for c in cases)
 sequence=(ROOT/'validation/fixtures/primer-design-template.txt').read_text().strip()
 rows+='D\tdesign\t'+sequence+'\n'
 output=subprocess.run([str(probe)],input=rows,text=True,capture_output=True,check=True).stdout
 actual={}; designs=[]
 for line in output.splitlines():
  f=line.split('\t')
  if f[0]=='T': actual[f[1]]=list(map(float,f[2:]))
  else: designs.append(f[3:])
 max_error=0.
 for identity,seq,dna,na,mg,dntp,dmso,formamide,temp,model,salt in cases:
  kw=dict(dna_conc=dna,mv_conc=na,dv_conc=mg,dntp_conc=dntp)
  tm=primer3.calc_tm(seq,**kw,dmso_conc=dmso,formamide_conc=formamide,tm_method=['breslauer','santalucia'][model],salt_corrections_method=['schildkraut','santalucia','owczarzy'][salt])
  hairpin=primer3.calc_hairpin(seq,**kw,temp_c=temp)
  homodimer=primer3.calc_homodimer(seq,**kw,temp_c=temp)
  expected=[tm,hairpin.tm,hairpin.dh,hairpin.ds,hairpin.dg,homodimer.tm,homodimer.dh,homodimer.ds,homodimer.dg]
  for found,wanted in zip(actual[identity],expected):
   max_error=max(max_error,abs(found-wanted))
   assert math.isclose(found,wanted,abs_tol=1e-7,rel_tol=1e-9),(identity,found,wanted)
 stored=json.loads((ROOT/'validation/fixtures/primer3-design-settings.json').read_text())
 assert hashlib.sha256(sequence.encode()).hexdigest()==stored['template_sha256']
 reference=primer3.bindings.design_primers(dict(SEQUENCE_TEMPLATE=sequence,SEQUENCE_TARGET=stored['target']),stored['settings'])
 assert len(designs)==reference['PRIMER_PAIR_NUM_RETURNED']==5
 for i,row in enumerate(designs):
  fs,fb,fe,rs,rb,re,ft,rt,size,penalty=row
  left=reference[f'PRIMER_LEFT_{i}']; right=reference[f'PRIMER_RIGHT_{i}']
  assert (fs,int(fb),int(fe))==(reference[f'PRIMER_LEFT_{i}_SEQUENCE'],left[0],sum(left))
  assert (rs,int(rb),int(re))==(reference[f'PRIMER_RIGHT_{i}_SEQUENCE'],right[0]-right[1]+1,right[0]+1)
  for found,key in [(ft,f'PRIMER_LEFT_{i}_TM'),(rt,f'PRIMER_RIGHT_{i}_TM'),(size,f'PRIMER_PAIR_{i}_PRODUCT_SIZE'),(penalty,f'PRIMER_PAIR_{i}_PENALTY')]:
   assert math.isclose(float(found),reference[key],abs_tol=1e-7),(i,key,found,reference[key])
 independent=MeltingTemp.Tm_NN(cases[0][1],nn_table=MeltingTemp.DNA_NN3,dnac1=25,dnac2=25,Na=50,Mg=1.5,dNTPs=.6,saltcorr=5)
 default_id=next(c[0] for c in cases if c[-2:]==(1,1) and c[2:9]==(50,50,1.5,.6,0,0,37) and c[1]==cases[0][1])
 assert abs(independent-actual[default_id][0])<0.02,(independent,actual[default_id][0])
 result=dict(result='pass',primer3_py=primer3.__version__,reaction_sequence_cases=len(cases),design_pairs_compared=len(designs),maximum_absolute_wrapper_error=max_error,independent_biopython_tm=independent,production_tm=actual[default_id][0],scope=__doc__.strip(),references=['https://primer3.org/manual.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC19045/','https://doi.org/10.1021/bi702363u'])
 path=ROOT/'validation/results/primer3-oracle.json'; path.parent.mkdir(exist_ok=True); path.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
if __name__=='__main__':main()
