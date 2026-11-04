#!/usr/bin/env bash
# Read-only offline validation; manifest paths are relative to untouched source/.
set -euo pipefail
item_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
python3 - "$item_dir" <<'PY'
from pathlib import Path
import hashlib, importlib.util, json, sys
item = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location('download', item/'script/download.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
rows = (item/'manifest.sha256').read_text().splitlines()
assert rows, 'empty manifest'
seen = set()
for line in rows:
    digest, name = line.split('  ', 1)
    assert name in module.ASSETS and name not in seen, 'unexpected or duplicate source path'
    content = (item/'source'/name).read_bytes()
    assert hashlib.sha256(content).hexdigest() == digest, 'hash mismatch: '+name
    module.validate(name, content)
    seen.add(name)
assert seen == set(module.ASSETS), 'missing original asset'
assert {p.name for p in (item/'source').iterdir()} == seen, 'unmanifested source asset'
log = json.loads((item/'author_log.json').read_text())
assert log['dataset_id'] == 'sirenko_2013'
assert {r['path'] for r in log['sources']} == seen
for row in log['sources']:
    content = (item/'source'/row['path']).read_bytes()
    assert row['sha256'] == hashlib.sha256(content).hexdigest()
    assert row['bytes'] == len(content) and row['url'] == module.ASSETS[row['path']]
print('PASS: Sirenko-2013 original assets, formats, hashes and provenance')
PY
