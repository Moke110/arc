#!/usr/bin/env python3
"""Re-materialize pinned, DEG-independent drug evidence without changing source/."""
import hashlib
import json
from pathlib import Path
import urllib.request


def main():
    item = Path(__file__).resolve().parent.parent
    files = json.loads((item / 'script/tox_evidence_sources.json').read_text())
    for name, spec in files.items():
        path = item / 'local/tox_evidence' / Path(name).name
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            with urllib.request.urlopen(spec['url']) as response:
                payload = response.read()
            if hashlib.sha256(payload).hexdigest() != spec['sha256']:
                raise SystemExit('independent evidence download changed: ' + name)
            path.write_bytes(payload)
        if hashlib.sha256(path.read_bytes()).hexdigest() != spec['sha256']:
            raise SystemExit('independent evidence snapshot changed: ' + name)
    print('DToxS independent evidence snapshots PASS')


if __name__ == '__main__':
    main()
