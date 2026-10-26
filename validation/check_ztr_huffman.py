#!/usr/bin/env python3
"""Generate STHUFF77 vectors with the public BSD Staden C encoder, then check
actual production Rust decoding/import. Includes inline/interlaced tables,
predefined sets1–3, shared HUFF chunks and multiple blocks. No reference binary
is retained. Input source/header revisions are SHA256-pinned before compilation.
"""
import argparse
import ctypes as C
import hashlib
import json
import random
import subprocess
import tempfile
import zlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DIRECTORY=ROOT/"validation/fixtures/trace/extended"
SOURCE_SHA="45bbfe4a8076c9f1210116679efef687ad5b9e8b9c191eb17fbfdcfad7f5ca36"
HEADER_SHA="98dcb6349a18633370c01d547b152caf00c81a4d49c0ef563a4e752332f58e95"
WRAPPER=r'''
#include <stdlib.h>
#include <string.h>
#include "io_lib/deflate_interlaced.h"
static char *take(block_t *blk,int *length) {
    *length=blk->byte+(blk->bit!=0);char *data=(char *)blk->data;block_destroy(blk,1);return data;
}
char *encode_test(unsigned char *raw,int len,int cset,int records,int *length) {
    huffman_codeset_t *codes=generate_code_set(cset,records,raw,len,1,MAX_CODE_LEN,0);
    block_t *blk=block_create(NULL,2);unsigned char head[2]={77,cset};store_bytes(blk,head,2);
    store_codes(blk,codes,1);
    if(cset!=0){blk->byte=2;memset(blk->data+2,0,blk->alloc-2);}
    huffman_multi_encode(blk,codes,cset,raw,len);
    if(cset==0||cset>=128)huffman_codeset_destroy(codes);
    return take(blk,length);
}
char *header_test(unsigned char *raw,int len,int cset,int records,int *length) {
    huffman_codeset_t *codes=generate_code_set(cset,records,raw,len,1,MAX_CODE_LEN,0);
    block_t *blk=block_create(NULL,2);unsigned char head[2]={0,cset};store_bytes(blk,head,2);
    store_codes(blk,codes,1);huffman_codeset_destroy(codes);return take(blk,length);
}
char *blocks_test(unsigned char *raw,int len,int split,int unused,int *length) {
    block_t *blk=block_create(NULL,2);unsigned char head[2]={77,0};store_bytes(blk,head,2);
    huffman_codeset_t *first=generate_code_set(0,1,raw,split,1,MAX_CODE_LEN,0);
    store_codes(blk,first,0);huffman_multi_encode(blk,first,0,raw,split);huffman_codeset_destroy(first);
    huffman_codeset_t *second=generate_code_set(0,3,raw+split,len-split,1,MAX_CODE_LEN,0);
    store_codes(blk,second,1);huffman_multi_encode(blk,second,0,raw+split,len-split);huffman_codeset_destroy(second);
    return take(blk,length);
}
void free_test(void *p){free(p);}
'''
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-source",type=Path,default=Path("/tmp/popgene-deflate-interlaced.c"))
    parser.add_argument("--reference-header",type=Path,default=Path("/tmp/popgene-deflate-interlaced.h"))
    args=parser.parse_args();source=args.reference_source.read_bytes();header=args.reference_header.read_bytes()
    assert hashlib.sha256(source).hexdigest()==SOURCE_SHA,"Unexpected C encoder revision"
    assert hashlib.sha256(header).hexdigest()==HEADER_SHA,"Unexpected C header revision"
    vectors=[];rng=random.Random(771026)
    with tempfile.TemporaryDirectory(prefix="popgene-sthuff-oracle-") as temporary:
        temporary=Path(temporary);(temporary/"io_lib").mkdir();(temporary/"io_lib/deflate_interlaced.h").write_bytes(header)
        (temporary/"deflate.c").write_bytes(source);(temporary/"wrapper.c").write_text(WRAPPER)
        library=temporary/"reference.dylib"
        subprocess.run(["cc","-shared","-fPIC","-O2","-I",str(temporary),str(temporary/"deflate.c"),str(temporary/"wrapper.c"),"-o",str(library)],check=True)
        codec=C.CDLL(str(library));codec.free_test.argtypes=[C.c_void_p]
        def call(name,raw,cset,records):
            function=getattr(codec,name);function.restype=C.c_void_p;length=C.c_int()
            data=function(C.c_char_p(raw),C.c_int(len(raw)),C.c_int(cset),C.c_int(records),C.byref(length));assert data
            output=C.string_at(data,length.value);codec.free_test(data);return output
        def record(name,raw,cset,records):
            table=call("header_test",raw,cset,records) if cset>=128 else None
            encoded=call("encode_test",raw,cset,records)
            vectors.append(dict(name=name,input=encoded.hex(),expected=raw.hex(),table=table.hex() if table else None))
        inputs=[b"\0",b"\0ACGTN"*100,b"\0English text, punctuation! 12345\n"*20,b"\0"+bytes(range(256))*3,b"\0"+bytes(rng.randrange(256) for _ in range(2048))]
        for cset in [0,1,2,3]:
            for i,raw in enumerate(inputs):record(f"set{cset}-{i}",raw,cset,1)
        for records in [2,3,4,8,16,64,256]:
            raw=b"\0"+bytes(rng.randrange(256) for _ in range(records*17+3))
            record(f"inline-interlaced{records}",raw,0,records)
            record(f"shared-interlaced{records}",raw,128,records)
        for split in [1,17,255]:
            raw=inputs[-1];vectors.append(dict(name=f"multiple-blocks-{split}",input=call("blocks_test",raw,split,0).hex(),expected=raw.hex(),table=None))
        bases=b"ACGTACGT";n=len(bases)
        channels=[[c*500+i*101 for i in range(n)] for c in range(4)]
        signal=b"\0\0"+b"".join(v.to_bytes(2,"big") for row in channels for v in row)
        payloads=[(b"BASE",b"\0"+bases,128,1),(b"SMP4",signal,129,4),(b"CNF1",b"\0"+bytes(range(21,29)),130,1)]
        def chunk(kind,data):return kind+b"\0"*4+len(data).to_bytes(4,"big")+data
        data=b"\xaeZTR\r\n\x1a\n\1\2"
        for kind,raw,cset,records in payloads:data+=chunk(kind,call("encode_test",raw,cset,records))
        data+=chunk(b"BPOS",call("encode_test",b"\0"*4+b"".join(i.to_bytes(4,"big") for i in range(n)),0,4))
        for _,raw,cset,records in payloads:
            table=call("header_test",raw,cset,records)
            if cset==128:table=call("encode_test",table,3,1)
            if cset==129:table=b"\2"+len(table).to_bytes(4,"little")+zlib.compress(table)
            data+=chunk(b"HUFF",table)
        DIRECTORY.mkdir(parents=True,exist_ok=True);(DIRECTORY/"shared-huffman.ztr").write_bytes(data)
    facts=dict(reference_sha256=SOURCE_SHA,reference_header_sha256=hashlib.sha256(header).hexdigest(),vectors=vectors,fixture=dict(name="shared-huffman.ztr",sha256=hashlib.sha256(data).hexdigest(),sequence=bases.decode(),qualities=list(range(21,29)),channels=channels))
    (DIRECTORY/"huffman-vectors.json").write_text(json.dumps(facts,indent=2)+"\n")
    subprocess.run([str(ROOT/"scripts/cargo-local"),"test","--locked","--offline","-p","popgene-formats","independent_staden_huffman","--lib"],cwd=ROOT,check=True)
    result=dict(status="pass",production_layer_cases=len(vectors),production_file_cases=1,reference_sha256=SOURCE_SHA,scope=__doc__.strip(),primary_references=["https://github.com/jkbonfield/io_lib/blob/master/io_lib/deflate_interlaced.c","https://www.rfc-editor.org/rfc/rfc1951"])
    (ROOT/"validation/results/ztr-huffman-oracle.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2))
if __name__=="__main__":main()
