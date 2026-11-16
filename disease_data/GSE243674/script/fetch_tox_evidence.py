#!/usr/bin/env python3
"""Re-materialize declared independent regulatory evidence; never touches source/."""
from pathlib import Path
import hashlib, urllib.request
FILES = [('ELLENCE_2025.pdf', 'https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/050778Orig1s029lbl.pdf', '76aaa8bec430e1355e5055c0a2243e9ab108227011fa35e1a309726839bd947c'), ('daunorubicin_DailyMed.html', 'https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=f4312095-3c21-296b-e053-2995a90a0d9c&type=display', 'b9e03e9833f242ea3e294bf75881a0a7b032e64e0b94377dd03cc794d249dbb9')]
out = Path(__file__).resolve().parent.parent / 'local/tox_evidence'
out.mkdir(parents=True, exist_ok=True)
for name, url, expected in FILES:
    path = out / name
    if not path.exists():
        with urllib.request.urlopen(url) as response:
            path.write_bytes(response.read())
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit('independent regulatory snapshot changed: ' + name)
