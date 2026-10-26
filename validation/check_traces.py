#!/usr/bin/env python3
"""Check production trace codecs against independent SCF and Staden-C ZTR facts.

The committed ZTR fixture was encoded by Staden io_lib ichebcomp in a temporary
C library, then packed and zlib-compressed independently in Python. Its source
hash and factual signal arrays are recorded in fixtures/trace/expected.json.
No vendor application code or executable is included.
"""
from pathlib import Path
import hashlib,json,subprocess,tempfile
from rust_probe import compile_probe
ROOT=Path(__file__).resolve().parents[1]
FIXTURES=ROOT/'validation/fixtures/trace'
RUST=r'''
use popgene_formats::{open_records,encode,import,SaveFormat};
use std::path::Path;
fn report(label:&str,document:&popgene_core::Document){
if let Some(trace)=&document.trace {let channels:Vec<_>=b"ACGT".iter().map(|base|&trace.channels[trace.channel_order.as_bytes().iter().position(|b|b==base).unwrap()]).collect();println!("{}\t{}\t{:?}\t{:?}\t{:?}",label,document.sequence,trace.qualities,trace.peaks,channels);}
else {println!("{}\t{}\t{:?}\t[]\t[]",label,document.sequence,document.sequence_qualities.as_ref().unwrap().iter().map(|q|q.expect("fixture qualities are known")).collect::<Vec<_>>());}}
fn main(){for filename in std::env::args().skip(1){let path=Path::new(&filename);let name=path.file_name().unwrap().to_str().unwrap();for (index,document) in open_records(path).unwrap().iter().enumerate(){report(&format!("{name}:{index}"),document);let formats=if document.trace.is_some(){vec![SaveFormat::Scf,SaveFormat::Ztr,SaveFormat::Fastq]}else{vec![SaveFormat::Fastq]};for format in formats{let data=encode(document,format).unwrap();let name=format!("roundtrip.{}",format.extension());let read=import(&name,&data).unwrap();report(&format!("{}:{index}:{}",path.file_name().unwrap().to_str().unwrap(),format.extension()),&read);}}}}
'''
def main():
    facts=json.loads((FIXTURES/'expected.json').read_text())
    for name,digest in facts['fixtures'].items():assert hashlib.sha256((FIXTURES/name).read_bytes()).hexdigest()==digest,name
    paths=sorted(p for p in FIXTURES.iterdir() if p.suffix in ('.scf','.ztr','.fastq'))
    with tempfile.TemporaryDirectory(prefix='popgene-trace-oracle-') as directory:
        directory=Path(directory);source=directory/'probe.rs';binary=directory/'probe';source.write_text(RUST)
        compile_probe(source,binary,['popgene-formats','popgene-core'])
        lines=subprocess.run([str(binary),*map(str,paths)],text=True,capture_output=True,check=True).stdout.splitlines()
    for line in lines:
        label,sequence,qualities,peaks,channels=line.split('\t');qualities=json.loads(qualities);peaks=json.loads(peaks);channels=json.loads(channels)
        if label.startswith('quality-reads.fastq'):
            index=int(label.split(':')[1]);assert sequence==['ACGT','TTTTGG'][index];assert qualities==[[0,20,40,93],[20]*6][index];assert not peaks and not channels
        else:
            assert sequence==facts['sequence'],label;assert qualities==facts['qualities'],label
            if not label.endswith(':fastq'):assert peaks==facts['peaks'],label;assert channels==facts['channels'],label
    result={'status':'pass','production_import_export_records':len(lines),'independent_formats':['SCF2','SCF3','ZTR1.2 integer-Chebyshev + signed16 packing + zlib','FASTQ Phred+33'],'staden_reference_codec_sha256':facts['reference_codec_sha256'],'exact_signals_per_trace':sum(map(len,facts['channels']))}
    output=ROOT/'validation/results/trace-oracle.json';output.parent.mkdir(exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
