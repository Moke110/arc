"""Archive the public Huang article/supplement; this is not the full screen."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import pymupdf

ITEM = Path(__file__).resolve().parent.parent
ROOT = next(p for p in ITEM.parents if (p / '.git').exists())
sys.path.insert(0, str(ROOT / 'script'))
from phenotype_archive import source_tree_summary, verify

NAME = 'mmc2 (1).pdf'
URL = 'https://www.axionbiosystems.com/sites/default/files/resources/mmc2%20%281%29.pdf'
ASSETS = {NAME: URL}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify or (ITEM / 'author_log.json').exists():
        verify(ITEM, 'huang_2022', None, ASSETS)
        return
    with tempfile.TemporaryDirectory() as temporary:
        payload = Path(temporary) / 'article.pdf'
        subprocess.run(['curl', '-sS', '-L', '--fail', '--max-time', '45',
                        '--max-filesize', '25000000', URL, '-o', str(payload)], check=True)
        content = payload.read_bytes()
        if not content.startswith(b'%PDF-') or len(content) > 25_000_000:
            raise ValueError('Unexpected PDF payload')
        with pymupdf.open(payload) as document:
            text = '\n'.join(page.get_text() for page in document)
        for marker in ['110643', 'Population-based', 'Table S3', 'Data and code availability']:
            if marker not in text:
                raise ValueError('Missing expected article/supplement content: ' + marker)
        source = ITEM / 'source'
        source.mkdir(exist_ok=True)
        dest = source / NAME
        if dest.exists() and dest.read_bytes() != content:
            raise ValueError('Refusing to replace differing original bytes')
        if not dest.exists():
            dest.write_bytes(content)
    digest = hashlib.sha256(content).hexdigest()
    (ITEM / 'manifest.sha256').write_text(f'{digest}  {NAME}\n')
    log = {
        'dataset_id': 'huang_2022',
        'upstream_identifier': '10.1016/j.celrep.2022.110643',
        'version': 'Public combined article and supplemental information from Axion',
        'retrieved_at': datetime.now(timezone.utc).isoformat(),
        'license': {'status': 'Article states CC BY-NC-ND 4.0; no full-dataset license established'},
        'validation': {'status': 'pdf_identity_and_hash_validated',
                       'scope': 'Article and public supplement only; full numerical screen not acquired'},
        'sources': [{'path': NAME, 'url': URL, 'media_type': 'application/pdf',
                     'bytes': len(content), 'sha256': digest}],
        'source_tree': source_tree_summary(ITEM),
        'local_additions': [],
    }
    (ITEM / 'author_log.json').write_text(json.dumps(log, indent=2) + '\n')
    print('Archived public article/supplement; full screen remains unavailable.')

if __name__ == '__main__':
    main()
