"""Fetch pinned upstream bytes without replacing an existing source file."""
import hashlib,json,urllib.request
from pathlib import Path
ITEM=Path(__file__).resolve().parents[1]
FILES=[{'name': 'GSE329393_raw_counts.tsv.gz', 'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329393/suppl/GSE329393_raw_counts.tsv.gz', 'sha256': '8dbee713d12059b279ed02cc2ec62af639c1f8baaae15eed121b710b1ea85f90', 'bytes': 591619}, {'name': 'GSE329393_series_matrix.txt.gz', 'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329393/matrix/GSE329393_series_matrix.txt.gz', 'sha256': '98570415473a718bf42b3eb002f51002bef0863b6ad2a27e0486233ce6a0db89', 'bytes': 4173}]
for r in FILES:
 p=ITEM/'source'/r['name'];p.parent.mkdir(exist_ok=True)
 if p.exists():
  if hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise SystemExit('existing source digest mismatch: '+str(p))
 else:
  b=urllib.request.urlopen(r['url'],timeout=120).read()
  if hashlib.sha256(b).hexdigest()!=r['sha256']:raise SystemExit('upstream digest changed')
  p.write_bytes(b)
