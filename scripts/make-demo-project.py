"""Assemble a disposable project from PopGene's original synthetic fixtures."""
from pathlib import Path
import shutil
root=Path(__file__).resolve().parents[1]
destination=root/'validation/projects/workspace-smoke'
files={'DNA Files/parent.dna':'annotated.dna','DNA Files/pcr-product.dna':'synthetic-pcr.dna','DNA Files/annotated.gb':'annotated.gb','Sequencing/synthetic.ab1':'synthetic-complete.ab1','Validation/truncated.dna':'truncated.dna'}
for target,source in files.items():
    path=destination/target
    path.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(root/'validation/fixtures'/source,path)
(destination/'README.txt').write_text('PopGene original synthetic project. No laboratory data. DNA and trace files should open in separate tabs. Validation/truncated.dna is intentionally malformed and must produce a per-file error without closing the other tabs. README.txt is an unsupported project attachment and must remain visible.\n')
print(destination)
