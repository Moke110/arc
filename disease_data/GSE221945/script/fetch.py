"""Fetch pinned upstream bytes without replacing an existing source file."""
import hashlib,json,urllib.request
from pathlib import Path
ITEM=Path(__file__).resolve().parents[1]
FILES=[{'name': 'GSE221945_Isogenic_vs_Brugada_min5_DESeq2.csv.gz', 'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221945/suppl/GSE221945_Isogenic_vs_Brugada_min5_DESeq2.csv.gz', 'sha256': 'a1f4ebf5b3dfa8435c9f79fe550fda36cc8b8013266a6574d035fedbb948fe7a', 'bytes': 1717914}, {'name': 'GSE221945_series_matrix.txt.gz', 'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221945/matrix/GSE221945_series_matrix.txt.gz', 'sha256': '7b45dce7adab6805fb0b73e828dc7435009057eb7c4e2d4ece50e60273b66b51', 'bytes': 3784}]
for r in FILES:
 p=ITEM/'source'/r['name'];p.parent.mkdir(exist_ok=True)
 if p.exists():
  if hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise SystemExit('existing source digest mismatch: '+str(p))
 else:
  b=urllib.request.urlopen(r['url'],timeout=120).read()
  if hashlib.sha256(b).hexdigest()!=r['sha256']:raise SystemExit('upstream digest changed')
  p.write_bytes(b)
