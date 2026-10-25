"""Sample the running local PopGene bundle; RSS is not peak/physical/GPU footprint."""
import argparse, hashlib, json, platform, subprocess, time
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('--label', required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
binary = root/'dist/PopGene.app/Contents/MacOS/PopGene'
bundle = binary.parents[2]
pids = subprocess.check_output(['pgrep', '-x', 'PopGene'], text=True).split()
matched=[]
for pid in pids:
    command=subprocess.check_output(['ps','-p',pid,'-o','comm='],text=True).strip()
    if command==str(binary): matched.append(pid)
if len(matched)!=1:
    raise SystemExit(f'Expected one running bundle, found {len(matched)}')
samples=[]
for i in range(5):
    rss,cpu=subprocess.check_output(['ps','-p',matched[0],'-o','rss=,pcpu='],text=True).split()
    samples.append({'rss_kib':int(rss),'cpu_percent':float(cpu)})
    if i<4: time.sleep(1)
result={'label':a.label,'observed_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'platform':platform.platform(),'machine':platform.machine(),
        'binary_bytes':binary.stat().st_size,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'bundle_file_bytes':sum(f.stat().st_size for f in bundle.rglob('*') if f.is_file()),
        'samples':samples,'method':'Five ps RSS samples, one second apart, after UI ready. Single application process; excludes separately accounted GPU/driver memory. Not peak RSS, physical footprint, cold-start timing, or a portable performance claim.'}
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
