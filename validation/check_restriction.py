#!/usr/bin/env python3
"""Compare production Rust with independent Biopython 1.88 regexes/offset data.
Run .venv/bin/python validation/check_restriction.py. No network is used.
"""
import json
from pathlib import Path
import random
import re
import subprocess
import Bio
from rust_probe import build_example
from Bio import Restriction, SeqIO
from Bio.Seq import Seq
ROOT = Path(__file__).resolve().parents[1]
IUPAC = dict(zip('ACGTRYSWKMBDHVN', ['A','C','G','T','AG','CT','CG','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT']))
PRIMARY = {
 'EcoP15I':'https://www.neb.com/en/products/r0646-ecop15i',
 'HpyUM037X':'https://www.neb.com/ja-jp/-/media/nebjp/catalog-and-literature-request/id3_jp_lr_restriction_enzyme_handbook.pdf',
 'catalogue':'https://github.com/biopython/biopython/blob/biopython-188/Bio/Restriction/Restriction_Dictionary.py',
}
def main():
 assert Bio.__version__ == '1.88', Bio.__version__
 probe = build_example('popgene-core', 'restriction_probe')
 def run(rows): return subprocess.run([str(probe)],input=rows,text=True,capture_output=True,check=True).stdout
 catalogue = {}
 for line in run('').splitlines():
  _,name,motif,known,top,bottom,second_top,second_bottom = line.split('\t')
  cuts = tuple(None if v == '-' else int(v) for v in (top,bottom,second_top,second_bottom))
  catalogue[name] = motif,known=='true',cuts
  if name == 'EcoP15I': assert (motif,known,cuts)==('CAGCAG','true',(31,33,None,None))
  else:
   e = getattr(Restriction,name)
   assert motif == e.site
   assert (known=='true') == (e.fst5 is not None and e.fst3 is not None)
   expected = (e.fst5 or 0,e.size+e.fst3 if e.fst3 is not None else 0,e.scd5,e.size+e.scd3 if e.scd3 is not None else None)
   assert cuts == expected,(name,cuts,expected)
 # Evaluate Bio's independent forward/reverse regex groups separately. A single
 # alternation hides the reverse orientation when both match at the same start.
 regexes = {}
 for name,(motif,_,_) in catalogue.items():
  if name in ('EcoP15I','HpyUM037X'):
   # EcoP15I absent from Bio; Hpy's Bio regex omits its documented second motif.
   pairs = []
   for alternative in motif.split('|'):
    for reverse,s in [(False,alternative),(True,str(Seq(alternative).reverse_complement()))]:
     pattern = ''.join('['+IUPAC[b]+']' for b in s)
     pairs.append((reverse,re.compile('(?=('+pattern+'))')))
   regexes[name] = pairs
  else:
   groups = re.findall(r'\(\?P<([^>]+)>([^)]+)\)',getattr(Restriction,name).compsite.pattern)
   regexes[name] = [(label.endswith('_as'),re.compile('(?=('+pattern+'))')) for label,pattern in groups]
   assert regexes[name],name
 rng = random.Random(823)
 dna = lambda length: ''.join(rng.choices('ACGT',k=length))
 cases = []
 fixture = str(SeqIO.read(ROOT/'validation/fixtures/annotated.gb','genbank').seq)
 for circular in (False,True): cases.append((f'annotated-{circular}',circular,'*',fixture))
 for kind,seq in [('random',dna(4000)),('lowercase',dna(1200).lower()),('ambiguous','N'.join(dna(19) for _ in range(60)))]:
  for circular in (False,True): cases.append((f'{kind}-{circular}',circular,'*',seq))
 for name,(motif,_,_) in catalogue.items():
  for ai,alternative in enumerate(motif.split('|')):
   for reverse in (False,True):
    concrete = ''.join(rng.choice(IUPAC[b]) for b in alternative)
    if reverse: concrete = str(Seq(concrete).reverse_complement())
    sequence = dna(80)+concrete+dna(80)
    prefix = f'{name}-{ai}-{reverse}'
    cases.append((prefix+'-linear',False,name,sequence))
    rotation = 80+len(concrete)//2
    cases.append((prefix+'-circular',True,name,sequence[rotation:]+sequence[:rotation]))
 for name in ['EcoRI','BcgI','SapI','DpnII','HpyUM037X']:
  concrete = ''.join(IUPAC[b][0] for b in catalogue[name][0].split('|')[0])
  cases.append((name+'-edge',False,name,concrete))
 counts,sites = {},{}
 rows = ''.join(f"{name}\t{'circular' if circle else 'linear'}\t{enzyme}\t{seq}\n" for name,circle,enzyme,seq in cases)
 for line in run(rows).splitlines():
  f = line.split('\t')
  if f[0]=='C': counts[f[1],f[2]]=int(f[3])
  elif f[0]=='S': sites.setdefault((f[1],f[2]),[]).append((int(f[3]),f[4]=='true',int(f[5]),*(None if v=='-' else int(v) for v in f[6:])))
 comparisons = 0
 for case,circular,selection,sequence in cases:
  n = len(sequence)
  for name in catalogue if selection=='*' else [selection]:
   motif,known,cuts = catalogue[name]
   scan = sequence.upper()+(sequence.upper()[:max(map(len,motif.split('|')))-1] if circular else '')
   expected = []
   for reverse,regex in regexes[name]:
    for hit in regex.finditer(scan):
     recognized = hit.group(1)
     if hit.start()>=n or len(recognized)>n or not set(recognized)<=set('ACGT'): continue
     start,length = hit.start(),len(recognized)
     def boundary(offset):
      if offset is None: return None
      cut = start+offset
      return cut%n if circular else cut if 0<=cut<=n else None
     offsets = cuts if known else (None,None,None,None)
     if reverse: offsets=tuple(None if v is None else length-v for v in [offsets[1],offsets[0],offsets[3],offsets[2]])
     expected.append((start,reverse,length,*(boundary(x) for x in offsets)))
   expected.sort(key=lambda s:(s[0],s[1]))
   assert counts[case,name]==len(expected),(case,name,counts[case,name],len(expected))
   assert sites.get((case,name),[])==expected[:512],(case,name,sites.get((case,name)),expected)
   comparisons += 1
 result = {
  'result':'pass','biopython_version':Bio.__version__,'catalogue_size':len(catalogue),
  'defined_cleavage':sum(known for _,known,_ in catalogue.values()),
  'two_cut_enzymes':sum(cuts[2] is not None for _,_,cuts in catalogue.values()),
  'cases':len(cases),'enzyme_case_comparisons':comparisons,
  'scope':'Actual Rust scanner versus independent Biopython regex matching and catalogue offsets: both cut pairs, unknown cuts, case, ambiguity, origin crossing and molecule edges.',
  'oracle_exceptions':['EcoP15I primary-source motif/offsets; absent from Biopython.','HpyUM037X primary-source second motif GTGGNAG; omitted by Biopython regex but present in site data.','Each regex orientation is tested independently; Biopython alternation can hide reverse matches when both orientations occur at one start.'],
  'primary_references':PRIMARY,
  'fixture_ecori':sites.get(('annotated-True','EcoRI')),'fixture_nmeaiii':sites.get(('annotated-True','NmeAIII')),
 }
 path = ROOT/'validation/results/restriction-oracle.json'
 path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__': main()
