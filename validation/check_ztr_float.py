#!/usr/bin/env python3
"""Validate legacy ZTR73 against the published disabled Staden C encoder.
Pins binary32 arithmetic, disables FMA contraction, fixes only its short-block
returned-length bug. Decoder bounds bugs are not copied. Historical encoders
using x87 excess precision/libm differences remain an explicitly warned risk.
"""
import ctypes as C,hashlib,json,random,subprocess,tempfile
from pathlib import Path
from check_ztr_extended import function,SHA256,URL
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'validation/fixtures/trace/floating73'
PREAMBLE='''
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#define ZTR_FORM_CHEB445 73
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define be_int2(x) ((uint16_t)__builtin_bswap16((uint16_t)(x)))
#define le_int2(x) ((uint16_t)(x))
#else
#define be_int2(x) ((uint16_t)(x))
#define le_int2(x) ((uint16_t)__builtin_bswap16((uint16_t)(x)))
#endif
'''
def main():
    source=Path('/tmp/popgene-iolib-compression.c').read_bytes();assert hashlib.sha256(source).hexdigest()==SHA256
    encoder=function(source.decode(),'cheb445comp').replace('*data_len = nwords*2;', '*data_len = (nwords+1)*2;')
    rng=random.Random(73026);vectors=[]
    with tempfile.TemporaryDirectory(prefix='popgene-floating-ztr-') as temporary:
        temporary=Path(temporary);cfile=temporary/'reference.c';library=temporary/'reference.dylib'
        cfile.write_text(PREAMBLE+encoder+'\nvoid release(void *data){free(data);}\n')
        subprocess.run(['cc','-shared','-fPIC','-O2','-ffp-contract=off',str(cfile),'-lm','-o',str(library)],check=True)
        codec=C.CDLL(str(library));codec.cheb445comp.restype=C.c_void_p;codec.release.argtypes=[C.c_void_p]
        def encode(raw):
            length=C.c_int();ptr=codec.cheb445comp(C.c_char_p(raw),C.c_int(len(raw)),C.byref(length));assert ptr
            result=C.string_at(ptr,length.value);codec.release(ptr);return result
        for index in range(160):
            n=index+1 if index<6 else rng.randrange(5,130)
            values=[rng.randrange(256) if index%2 else 0]+[rng.randrange(65536) if index%2 else (i*i*7)%65536 for i in range(n-1)]
            raw=b''.join(v.to_bytes(2,'big') for v in values);data=encode(raw)
            vectors.append(dict(name=f'C-binary32-{index}-words{n}',input=data.hex(),expected=raw.hex()))
        channels=[[0,10,100,200,300,400],[1,2,3,4,5,6],[600,500,400,300,200,100],[20,40,60,80,100,120]]
        def chunk(tag,data):return tag+b'\0'*4+len(data).to_bytes(4,'big')+data
        ztr=b'\xaeZTR\r\n\x1a\n\1\2'
        ztr+=chunk(b'SMP4',encode(b'\0\0'+b''.join(v.to_bytes(2,'big') for channel in channels for v in channel)))
        ztr+=chunk(b'BASE',b'\0ACGT')+chunk(b'BPOS',b'\0'*4+b''.join(i.to_bytes(4,'big') for i in [1,2,3,4]))+chunk(b'CNF1',b'\0'+bytes([30]*4))
        DEST.mkdir(parents=True,exist_ok=True);(DEST/'legacy.ztr').write_bytes(ztr)
    facts={'source':URL,'source_sha256':SHA256,'semantics':'Little-endian reference first word; IEEE754 binary32 intermediates, double outer expression, truncation-to-int, no FMA contraction; encoder short-block returned length fixed to include header','vectors':vectors}
    (DEST/'vectors.json').write_text(json.dumps(facts,indent=2)+'\n')
    subprocess.run([str(ROOT/'scripts/cargo-local'),'test','--locked','--offline','-p','popgene-formats','floating_filter::tests'],cwd=ROOT,check=True)
    result={'status':'pass','independent_c_vectors':len(vectors),'full_ztr_imports':1,'source_sha256':SHA256,'semantics':facts['semantics']}
    (ROOT/'validation/results').mkdir(exist_ok=True);(ROOT/'validation/results/ztr-floating-oracle.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
