"""Original synthetic ABIF parameter fixture; contains no experimental data."""
from pathlib import Path
import struct
import hashlib
import json
from generate_fixtures import abif

ROOT = Path(__file__).resolve().parent
source = abif()
count, offset = struct.unpack_from('>I', source, 18)[0], struct.unpack_from('>I', source, 26)[0]
tags = []
for i in range(count):
    pos = offset + 28*i
    name, num, typ, size, elements, length = struct.unpack_from('>4sIHHII', source, pos)
    start = pos+20 if length <= 4 else struct.unpack_from('>I', source, pos+20)[0]
    tags.append((name,num,typ,size,source[start:start+length]))
def pstr(s):
    b=s.encode('ascii'); return bytes([len(b)])+b
tags += [
    (b'SMPL',1,18,1,pstr('PopGene synthetic sample')),
    (b'MCHN',1,18,1,pstr('SYNTHETIC-0001')),
    (b'MODL',1,2,1,b'3730'),
    (b'TUBE',1,18,1,pstr('A01')),
    (b'LANE',1,4,2,struct.pack('>h',1)),
    (b'RunN',1,19,1,b'Synthetic parity run\0'),
    (b'RUND',1,10,4,struct.pack('>HBB',2026,10,22)),
    (b'RUNT',1,11,4,bytes([9,30,15,0])),
    (b'CMNT',1,18,1,pstr('Original synthetic metadata, not an instrument run')),
    (b'DySN',1,18,1,pstr('Synthetic dye set')),
    (b'SPAC',1,7,4,struct.pack('>f',3.0)),
    (b'S/N%',1,4,2,struct.pack('>4h',10,20,30,40)),
    (b'SCAN',1,5,4,struct.pack('>i',16)),
    (b'SVER',1,18,1,pstr('Synthetic version 1.0')),
]
entries=[]; payload=bytearray(); payload_start=128+28*len(tags)
for name,num,typ,size,data in tags:
    field=data.ljust(4,b'\0') if len(data)<=4 else struct.pack('>I',payload_start+len(payload))
    entries.append(struct.pack('>4sIHHII',name,num,typ,size,len(data)//size,len(data))+field+struct.pack('>I',0))
    if len(data)>4: payload.extend(data)
root=struct.pack('>4sIHHIIII',b'tdir',1,1023,28,len(tags),28*len(tags),128,0)
result=(b'ABIF'+struct.pack('>H',101)+root).ljust(128,b'\0')+b''.join(entries)+payload
path=ROOT/'fixtures/regression/trace-parameters.ab1'
path.parent.mkdir(parents=True,exist_ok=True)
if path.exists() and path.read_bytes()!=result: raise SystemExit('Refusing to overwrite different fixture')
path.write_bytes(result)
print(json.dumps({'path':str(path),'sha256':hashlib.sha256(result).hexdigest(),'provenance':'Original synthetic calls, signals and ABIF tags; never an experimental read'},indent=2))
