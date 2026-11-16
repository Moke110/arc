"""Fetch pinned upstream bytes without replacing an existing source file."""
import hashlib,json,urllib.request
from pathlib import Path
ITEM=Path(__file__).resolve().parents[1]
FILES=[{'name': 'GSE248187_DSP_EHTs_Differential_Expression_Table.txt.gz', 'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248187/suppl/GSE248187_DSP_EHTs_Differential_Expression_Table.txt.gz', 'sha256': '46a79a7d684e539e3cd751eae81398e75d6d0a7b176a8077b29039f70dcf9bda', 'bytes': 1521590}, {'name': 'GSE248187.xml', 'url': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11213508/fullTextXML', 'sha256': '3aac314799acf90648dcf2bb06d68aaf2b079da3a47e2c5f776d2ee4d08f6e7d', 'bytes': 175730}]
for r in FILES:
 p=ITEM/'source'/r['name'];p.parent.mkdir(exist_ok=True)
 if p.exists():
  if hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise SystemExit('existing source digest mismatch: '+str(p))
 else:
  b=urllib.request.urlopen(r['url'],timeout=120).read()
  if hashlib.sha256(b).hexdigest()!=r['sha256']:raise SystemExit('upstream digest changed')
  p.write_bytes(b)
