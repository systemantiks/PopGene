#!/usr/bin/env python3
"""Generate ZTR vectors with the independent BSD Staden C reference functions,
then test the actual Rust decoder and complete-file importer. No application
binaries or reference shared library are retained. Supply compression.c with
--reference-source; the pinned public-source hash is checked before compilation.
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

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "validation/fixtures/trace/extended"
SHA256 = "fc677c17b29be0029caa95444829b99c2eace616f36b25e28af218026b560b52"
URL = "https://raw.githubusercontent.com/jkbonfield/io_lib/master/io_lib/compression.c"
PREAMBLE = r"""
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#define xmalloc malloc
#define xfree free
#define ZTR_FORM_XRLE 3
#define ZTR_FORM_XRLE2 4
#define ZTR_FORM_LOG2 75
#define ZTR_FORM_QSHIFT 79
#define ZTR_FORM_TSHIFT 80
#define ZTR_TYPE_BASE 0
typedef struct { char *data; int dlength; } ztr_chunk_t;
typedef struct { ztr_chunk_t base; } ztr_t;
ztr_chunk_t **ztr_find_chunks(ztr_t *ztr, int unused, int *count) {
    ztr_chunk_t **result=malloc(sizeof(*result)); *count=1; result[0]=&ztr->base; return result;
}
int uncompress_chunk(ztr_t *ztr, ztr_chunk_t *chunk) { return 0; }
"""
WRAPPERS = r"""
char *test_tshift(char *data,int length,char *base,int base_length,int *out_length) {
    ztr_t ztr={{base,base_length}};return tshift(&ztr,data,length,out_length);
}
char *test_untshift(char *data,int length,char *base,int base_length,int *out_length) {
    ztr_t ztr={{base,base_length}};return untshift(&ztr,data,length,out_length);
}
void test_free(void *data) { free(data); }
"""
def function(source,name):
    start=source.index("char *"+name+"(")
    pos=source.index("{",start);depth=1;end=pos+1
    while depth:
        depth+=(source[end]=="{")-(source[end]=="}");end+=1
    return source[start:end]

def main():
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument("--reference-source",type=Path,default=Path("/tmp/popgene-iolib-compression.c"))
    source=args.parse_args().reference_source.read_bytes()
    assert hashlib.sha256(source).hexdigest()==SHA256,"Unexpected reference source revision"
    vectors=[];rng=random.Random(261025)
    with tempfile.TemporaryDirectory(prefix="popgene-ztr-reference-") as temporary:
        temporary=Path(temporary);cfile=temporary/"codec.c";library=temporary/"codec.dylib"
        names=["xrle","unxrle","xrle2","unxrle2","log2_data","unlog2_data","qshift","unqshift","tshift","untshift"]
        cfile.write_text(PREAMBLE+"\n".join(function(source.decode(),n) for n in names)+WRAPPERS)
        # The reference XRLE decoder compares char payload bytes with an
        # unsigned guard. Use unsigned char explicitly so guards128..255 have
        # the byte semantics prescribed by the format on every C platform.
        subprocess.run(["cc","-shared","-fPIC","-funsigned-char","-O2",str(cfile),"-lm","-o",str(library)],check=True)
        codec=C.CDLL(str(library));codec.test_free.argtypes=[C.c_void_p]
        def call(name,data,*parameters):
            f=getattr(codec,name);f.restype=C.c_void_p
            args=[C.c_char_p(data),C.c_int(len(data))]
            for p in parameters:args.extend([C.c_char_p(p),C.c_int(len(p))] if isinstance(p,bytes) else [C.c_int(p)])
            length=C.c_int();pointer=f(*args,C.byref(length));assert pointer
            result=C.string_at(pointer,length.value);codec.test_free(pointer);return result
        def record(name,data,expected,bases=None):
            vectors.append(dict(name=name,input=data.hex(),expected=expected.hex(),bases=None if bases is None else bases.hex(),lossy=data[0]==75))
        for size in [1,2,4,8]:
            word=bytes(rng.randrange(256) for _ in range(size))
            raw=b"\0"+word*300+bytes(range(256))+word*2
            for guard in [0,255]:
                encoded=call("xrle",raw,guard,size)
                assert call("unxrle",encoded)==raw
                record(f"XRLE-size{size}-guard{guard}",encoded,raw)
            raw=(b"\0"*size)+word*300+bytes(rng.randrange(256) for _ in range(size*200))+word*2
            encoded=call("xrle2",raw,size)
            assert call("unxrle2",encoded)==raw
            record(f"XRLE2-size{size}",encoded,raw)
        # Test every representable LOG2 code, including 0 and the 65535 endpoint.
        encoded=bytes([75,0])+b"".join(v.to_bytes(2,"big") for v in range(161))
        record("LOG2-all-codes",encoded,call("unlog2_data",encoded))
        for n in [0,1,4,257]:
            raw=b"\0"+bytes(rng.randrange(256) for _ in range(n*4));encoded=call("qshift",raw)
            assert call("unqshift",encoded)==raw
            record(f"QSHIFT-{n}",encoded,raw)
        for bases in [b"ACGT",b"AcgTNn",b"TCGAGCAT"*20]:
            raw=b"\0\0"+bytes(rng.randrange(256) for _ in range(len(bases)*8))
            encoded=call("test_tshift",raw,b"\0"+bases)
            assert call("test_untshift",encoded,b"\0"+bases)==raw
            record("TSHIFT-"+bases[:8].decode()+f"-{len(bases)}",encoded,raw,bases)
        def chunk(kind,data):return kind+b"\0"*4+len(data).to_bytes(4,"big")+data
        def zlib_layer(data):return b"\2"+len(data).to_bytes(4,"little")+zlib.compress(data)
        bases=b"ACGTACGT";n=len(bases)
        raw=b"\0\0"+b"".join((i*173+channel*510).to_bytes(2,"big") for channel in range(4) for i in range(n))
        quality=b"\0"+bytes(range(21,21+n))+bytes(range(3*n))
        fixtures=[]
        for logarithmic in [False,True]:
            signal=call("test_tshift",raw,b"\0"+bases)
            expected=raw
            if logarithmic:
                encoded_log=call("log2_data",raw)
                expected=call("unlog2_data",encoded_log)
                signal=encoded_log
            else:signal=call("xrle2",signal,4)
            data=b"\xaeZTR\r\n\x1a\n\1\2"+chunk(b"SMP4",zlib_layer(signal))+chunk(b"CNF4",call("qshift",quality))+chunk(b"BPOS",b"\0"*4+b"".join(i.to_bytes(4,"big") for i in range(n)))+chunk(b"BASE",call("xrle",b"\0"+bases,255,2))
            filename="extended-log2.ztr" if logarithmic else "extended-shifts.ztr"
            DIRECTORY.mkdir(parents=True,exist_ok=True);(DIRECTORY/filename).write_bytes(data)
            samples=[int.from_bytes(expected[i:i+2],"big") for i in range(2,len(expected),2)]
            fixtures.append(dict(name=filename,sha256=hashlib.sha256(data).hexdigest(),sequence=bases.decode(),qualities=list(quality[1:n+1]),channels=[samples[c*n:(c+1)*n] for c in range(4)],lossy=logarithmic))
    facts=dict(reference_source=URL,reference_sha256=SHA256,reference_c_flags=["-funsigned-char"],vectors=vectors,fixtures=fixtures)
    (DIRECTORY/"vectors.json").write_text(json.dumps(facts,indent=2)+"\n")
    subprocess.run([str(ROOT/"scripts/cargo-local"),"test","--locked","--offline","-p","popgene-formats","independent_staden_extended","--lib"],cwd=ROOT,check=True)
    result=dict(status="pass",independent_reference="Staden io_lib C compression.c",reference_sha256=SHA256,methods=[3,4,75,79,80],production_layer_cases=len(vectors),production_file_cases=len(fixtures),log2="All codes0..160 checked; quantization explicitly reported")
    (ROOT/"validation/results/ztr-extended-oracle.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2))

if __name__=="__main__":main()
