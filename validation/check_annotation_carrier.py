"""Independent reader for the proposed public v1 comment convention.

Usage: python3 validation/check_annotation_carrier.py FILE [FILE ...]

Generate examples by setting POPGENE_LINEAGE_QA_DIR to a new empty directory
when running the formats lineage tests, then pass those files to this checker.
Uses only Python's standard library for container/text framing and hashes.
This checks exact published record shapes, checksum integrity, raw annotation
record binding and limited native annotation associations. It does NOT prove
context_sha256, restoration semantics, authenticity, or absent ancestor DNA.
The production core and independent construction oracle cover those separate
semantic checks. No file is modified.
"""
from pathlib import Path
import sys,json,hashlib,re,uuid,xml.etree.ElementTree as ET
START='[PopGene lineage v'; END='[PopGene lineage end]'; MAX=4*1024*1024
REGISTRY_KEYS={'version','instance_id','tracking_limited','feature_epochs','origins','pieces'}
ORIGIN_KEYS={'id','source_id','occurrence_id','top','bottom','feature','references','extent','molecule','circular','source_digest','sequence_digest'}
PIECE_KEYS={'mapping_digest','origin_id','copy_id','feature_id','axis','runs','baseline','references_digest','sequence_digest'}
RUN_KEYS={'leaf','offset','length','target_start','reverse','strands'}
FEATURE_KEYS={'id','label','kind','parts','reverse','qualifiers','edited'}
REFERENCE_KEYS={'title','authors','consortium','journal','pubmed','remark','scope','extra'}
def exact_keys(value,keys):
 assert isinstance(value,dict), 'expected a JSON object'
 assert set(value)==keys, ('unexpected fields',set(value)^keys)

def unsigned(value, bits=64):
 assert type(value) is int and 0 <= value < (1 << bits), 'expected an unsigned integer'

def boolean(value):
 assert type(value) is bool, 'expected a JSON boolean'

def string(value):
 assert isinstance(value,str), 'expected a JSON string'

def digest(value):
 string(value)
 assert re.fullmatch('[0-9a-f]{64}',value), 'expected a SHA-256 digest'

def identifier(value):
 string(value)
 uuid.UUID(value)

def array(value, limit=None):
 assert isinstance(value,list), 'expected a JSON array'
 assert limit is None or len(value)<=limit, 'array exceeds limit'

def segment_shape(segment):
 exact_keys(segment,{'start','end'})
 unsigned(segment['start']);unsigned(segment['end'])

def location_shape(location, depth=0):
 assert depth<=64, 'location nesting exceeds64'
 assert isinstance(location,dict) and len(location)==1, 'expected a tagged location'
 tag,value=next(iter(location.items()))
 if tag=='Range':
  exact_keys(value,{'start','end','before','after'})
  unsigned(value['start']);unsigned(value['end'])
  boolean(value['before']);boolean(value['after'])
 elif tag=='Between':
  exact_keys(value,{'left','right'})
  unsigned(value['left']);unsigned(value['right'])
 elif tag=='Complement':
  location_shape(value,depth+1)
 elif tag in ('Join','Order','Bond','OneOf'):
  array(value,100_000)
  for child in value:location_shape(child,depth+1)
 elif tag=='Remote':
  exact_keys(value,{'accession','location'});string(value['accession'])
  if value['location'] is not None:location_shape(value['location'],depth+1)
 elif tag=='Gap':
  exact_keys(value,{'length','unknown100'});boolean(value['unknown100'])
  if value['length'] is not None:unsigned(value['length'])
 else:
  raise AssertionError(('unknown location variant',tag))

def feature_shape(feature):
 assert isinstance(feature,dict), 'expected a feature object'
 assert set(feature) in (FEATURE_KEYS,FEATURE_KEYS|{'location'})
 unsigned(feature['id']);string(feature['label']);string(feature['kind'])
 boolean(feature['reverse']);boolean(feature['edited'])
 array(feature['parts'])
 for segment in feature['parts']:segment_shape(segment)
 # None is omitted by the canonical serializer, not an explicit null field.
 if 'location' in feature:location_shape(feature['location'])
 array(feature['qualifiers'])
 for qualifier in feature['qualifiers']:
  array(qualifier);assert len(qualifier)==2
  assert isinstance(qualifier[0],str) and (qualifier[1] is None or isinstance(qualifier[1],str))

def reference_shape(reference):
 exact_keys(reference,REFERENCE_KEYS);string(reference['title'])
 for key in REFERENCE_KEYS-{'title','extra'}:
  if reference[key] is not None:string(reference[key])
 array(reference['extra'])
 for pair in reference['extra']:
  array(pair);assert len(pair)==2
  string(pair[0]);string(pair[1])

def registry_shape(registry):
 exact_keys(registry,REGISTRY_KEYS)
 unsigned(registry['version'],32);assert registry['version']==1
 identifier(registry['instance_id']);boolean(registry['tracking_limited'])
 array(registry['feature_epochs'],8192)
 for pair in registry['feature_epochs']:
  array(pair);assert len(pair)==2
  unsigned(pair[0]);identifier(pair[1])
 array(registry['origins'],2048);array(registry['pieces'],8192)
 for origin in registry['origins']:
  exact_keys(origin,ORIGIN_KEYS)
  for key in ('id','occurrence_id','source_digest','sequence_digest'):digest(origin[key])
  identifier(origin['source_id']);unsigned(origin['extent'])
  string(origin['molecule']);assert origin['molecule'] in ('Dna','Rna','Protein')
  boolean(origin['circular'])
  segment_shape(origin['top']);segment_shape(origin['bottom'])
  feature_shape(origin['feature']);array(origin['references'])
  for reference in origin['references']:reference_shape(reference)
 runs=0
 for piece in registry['pieces']:
  exact_keys(piece,PIECE_KEYS)
  for key in ('mapping_digest','origin_id','references_digest','sequence_digest'):digest(piece[key])
  identifier(piece['copy_id']);unsigned(piece['feature_id'])
  string(piece['axis']);assert piece['axis'] in ('Top','Duplex')
  feature_shape(piece['baseline']);array(piece['runs'])
  runs+=len(piece['runs']);assert runs<=65536
  for run in piece['runs']:
   exact_keys(run,RUN_KEYS);unsigned(run['leaf'],32)
   for key in ('offset','length','target_start'):unsigned(run[key])
   unsigned(run['strands'],8);boolean(run['reverse'])
   assert run['length']>0 and run['strands'] in (1,2,3)

def depth_check(text):
 depth=0;quoted=False;escaped=False
 for c in text:
  if quoted:
   if escaped:escaped=False
   elif c=='\\':escaped=True
   elif c=='"':quoted=False
  elif c=='"':quoted=True
  elif c in '[{':
   depth+=1;assert depth<=96, 'JSON depth exceeds96'
  elif c in ']}':
   depth-=1;assert depth>=0
 assert depth==0 and not quoted

sha=lambda b:hashlib.sha256(b).hexdigest()
def read(path):
 data=path.read_bytes(); assert len(data)<=128*1024*1024
 record_hash=hashlib.sha256(); comments=[]; native_nodes=[]; sequence=None
 if data[:1]==b'\x09':
  format='native'; offset=0; packets=[]
  while offset<len(data):
   assert offset+5<=len(data)
   kind=data[offset];length=int.from_bytes(data[offset+1:offset+5],'big');offset+=5
   assert offset+length<=len(data)
   body=data[offset:offset+length];offset+=length;packets.append((kind,body));assert len(packets)<=10_000
   if kind in (5,10):
    record_hash.update(bytes([kind]));record_hash.update(length.to_bytes(8,'big'));record_hash.update(body)
    root=ET.fromstring(body)
    for node in root:
     if node.tag=='Feature':native_nodes.append(node)
     elif node.tag=='Primer':
      seen=set()
      for site in node.findall('BindingSite'):
       key=(site.get('location'),site.get('boundStrand'))
       if key not in seen:native_nodes.append((node,site));seen.add(key)
      if not seen:native_nodes.append((node,None))
   if kind==6:
    comments.extend(node.text or '' for node in ET.fromstring(body).findall('Comments'))
   if kind in (0,21,32):sequence=body[1:].decode().upper()
 else:
  text=data.decode(); format='genbank' if text.startswith('LOCUS') else 'embl'; in_features=False;in_comment=False
  for line in text.splitlines(keepends=True):
   if format=='genbank':
    if line.startswith('FEATURES'):in_features=True
    if line.startswith(('ORIGIN','CONTIG')):in_features=False
    if in_features:record_hash.update(line.encode())
    if line.startswith('COMMENT'):
     comments.append(line[12:].strip());in_comment=True
    elif in_comment and line.startswith(' '*12):comments[-1]+=' '+line[12:].strip()
    else:in_comment=False
   else:
    if line.startswith('FT'):record_hash.update(line.encode())
    if line.startswith('CC   '):
     if in_comment:comments[-1]+=' '+line[5:].strip()
     else:comments.append(line[5:].strip())
     in_comment=True
    else:in_comment=False
 joined='\n'.join(comments);assert joined.count(START)==joined.count(END)==1
 begin=joined.index(START);after=joined.index(' begin]',begin);assert joined[begin+len(START):after]=='1'
 end=joined.index(END,after);physical=joined[after+len(' begin]'):end]
 assert not any(c.isspace() and c not in ' \t\r\n\v\f' for c in physical)
 compact=''.join(c for c in physical if c not in ' \t\r\n\v\f')
 match=re.fullmatch(r'bytes=([0-9]+);sha256=([0-9a-f]{64});json=(.*)',compact,re.S);assert match
 payload=match[3].encode();assert len(payload)==int(match[1])<=MAX and sha(payload)==match[2]
 depth_check(match[3])
 value=json.loads(payload);exact_keys(value,{'format','records_sha256','context_sha256','registry'})
 digest(value['records_sha256']);digest(value['context_sha256'])
 assert value['format']==format and value['records_sha256']==record_hash.hexdigest()
 registry=value['registry'];registry_shape(registry)
 origins={origin['id']:origin for origin in registry['origins']};assert len(origins)==len(registry['origins'])
 for origin in origins.values():
  exact_keys(origin,ORIGIN_KEYS);feature_shape(origin['feature']);uuid.UUID(origin['source_id'])
  for key in ('id','occurrence_id','source_digest','sequence_digest'):assert re.fullmatch('[0-9a-f]{64}',origin[key])
 seen=set();runs=0
 for piece in registry['pieces']:
  exact_keys(piece,PIECE_KEYS);feature_shape(piece['baseline'])
  for run in piece['runs']:
   exact_keys(run,RUN_KEYS);assert run['length']>0 and run['strands'] in (1,2,3)
  key=(piece['axis'],piece['feature_id']);assert key not in seen;seen.add(key)
  uuid.UUID(piece['copy_id']);assert piece['origin_id'] in origins
  baseline=piece['baseline'];assert baseline['id']==piece['feature_id']
  runs+=len(piece['runs'])
  if format=='native':
   assert 1<=piece['feature_id']<=len(native_nodes), 'native annotation index is outside emitted records'
   node=native_nodes[piece['feature_id']-1]
   if isinstance(node,tuple):
    primer,site=node;assert primer.get('name')==baseline['label']
   else:assert node.get('name')==baseline['label'] and node.get('type')==baseline['kind']
 assert runs<=65536
 return {'path':str(path),'sha256':sha(data),'format':format,'payload_bytes':len(payload),'payload_sha256':sha(payload),'records_sha256':record_hash.hexdigest(),'instance_id':registry['instance_id'],'origins':len(origins),'pieces':len(seen),'runs':runs,'integrity_and_record_binding':'pass','authenticity_claimed':False,'exact_registry_origin_piece_run_keys':True,'recursive_typed_registry_shape_verified':True,'ancestor_document_or_sequence_fields_present':False,'semantic_context_digest_verified':False,'restoration_semantics_verified':False}

def main():
 if len(sys.argv)<2:
  raise SystemExit('Usage: python3 validation/check_annotation_carrier.py FILE [FILE ...]')
 results=[read(Path(name)) for name in sys.argv[1:]]
 print(json.dumps({'status':'pass','independent_parser':'Python standard library','files':results},indent=2))

if __name__=='__main__':main()
