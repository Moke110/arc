#!/usr/bin/env python3
"""Acquire pinned original Sirenko files; never rewrite a different source version."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import tempfile
import subprocess

ITEM = Path(__file__).resolve().parent.parent
ASSETS = {
    'publisher-mmc3.xls': 'https://ars.els-cdn.com/content/image/1-s2.0-S0041008X13004171-mmc3.xls',
    'article.xml': 'https://pmc-oa-opendata.s3.amazonaws.com/PMC3900303.1/PMC3900303.1.xml',
    'metadata.json': 'https://pmc-oa-opendata.s3.amazonaws.com/PMC3900303.1/PMC3900303.1.json',
}

def validate(name, content):
    if name.endswith('.xls') and not content.startswith(bytes.fromhex('d0cf11e0a1b11ae1')):
        raise ValueError('not a real OLE XLS: ' + name)
    if name.endswith('.xml'):
        import xml.etree.ElementTree as ET
        if ET.fromstring(content).tag != 'article':
            raise ValueError('not article XML')
    if name.endswith('.json'):
        if json.loads(content)['pmcid'] != 'PMC3900303':
            raise ValueError('wrong article metadata')

def main():
    manifest = ITEM / 'manifest.sha256'
    expected = {}
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            digest, name = line.split('  ', 1)
            expected[name] = digest
    source = ITEM / 'source'
    source.mkdir(exist_ok=True)
    records = []
    for name, url in ASSETS.items():
        dest = source / name
        if dest.exists():
            content = dest.read_bytes()
            media = {'publisher-mmc3.xls':'application/excel', 'article.xml':'application/xml', 'metadata.json':'application/json'}[name]
        else:
            with tempfile.TemporaryDirectory() as tmp:
                payload, headers = Path(tmp)/'payload', Path(tmp)/'headers'
                subprocess.run(['curl','-sS','-L','--fail','--max-time','30','-D',str(headers),url,'-o',str(payload)], check=True)
                content = payload.read_bytes()
                media = next((line.split(':',1)[1].strip() for line in reversed(headers.read_text().splitlines()) if line.lower().startswith('content-type:')), 'unknown')
        if len(content) > 10_000_000:
            raise ValueError('unexpected asset size')
        validate(name, content)
        digest = hashlib.sha256(content).hexdigest()
        if name in expected and digest != expected[name]:
            raise ValueError('pinned source hash differs: ' + name)
        if not dest.exists():
            fd, temp = tempfile.mkstemp(dir=source)
            with os.fdopen(fd, 'wb') as handle:
                handle.write(content)
            Path(temp).replace(dest)
        records.append({'path':name,'url':url,'media_type':media,'bytes':len(content),'sha256':digest})
    if not manifest.exists():
        manifest.write_text(''.join(f"{r['sha256']}  {r['path']}\n" for r in sorted(records,key=lambda r:r['path'])))
    log = ITEM / 'author_log.json'
    if not log.exists():
        log.write_text(json.dumps({'dataset_id':'sirenko_2013','upstream_identifier':'10.1016/j.taap.2013.09.017','version':'PMC3900303.1; publisher Table 3 bytes pinned by manifest','retrieved_at':datetime.now(timezone.utc).isoformat(),'sources':records,'license':{'evidence':'article.xml copyright and metadata.json license_code=TDM; is_pmc_openaccess=false','status':'rights reserved; public accessibility does not establish redistribution permission'},'validation':{'status':'validated_original_formats_and_sha256','table':'Supplemental Table 3, 24hr_beat_count'},'local_additions':[]},indent=2)+'\n')
    print('Sirenko-2013 acquired: '+str(len(records))+' pinned original assets')

if __name__ == '__main__':
    main()
