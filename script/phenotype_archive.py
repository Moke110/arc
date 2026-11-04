"""Shared bounded acquisition and offline validation for phenotype originals."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zipfile
import io


def validate(name, content, pmcid):
    """Reject download/error documents before originals enter the archive."""
    if name.endswith('.xls') and not content.startswith(bytes.fromhex('d0cf11e0a1b11ae1')):
        raise ValueError('not a real OLE XLS: ' + name)
    if name.endswith('.xlsx'):
        with zipfile.ZipFile(io.BytesIO(content)) as archive:
            if 'xl/workbook.xml' not in archive.namelist() or archive.testzip():
                raise ValueError('not an intact XLSX: ' + name)
    if name.endswith('.pdf') and not content.startswith(b'%PDF-'):
        raise ValueError('not PDF: ' + name)
    if name.endswith('.jpg') and not content.startswith(b'\xff\xd8\xff'):
        raise ValueError('not JPEG: ' + name)
    if name.endswith('.xml') and ET.fromstring(content).tag != 'article':
        raise ValueError('not article XML: ' + name)
    if name.endswith('.json') and json.loads(content)['pmcid'] != pmcid:
        raise ValueError('wrong article metadata')


def verify(item, dataset_id, pmcid, assets):
    """Check source-only inventory, bytes, formats and frozen provenance offline."""
    rows = (item/'manifest.sha256').read_text().splitlines()
    if not rows:
        raise ValueError('empty manifest')
    seen = set()
    for line in rows:
        digest, name = line.split('  ', 1)
        if name not in assets or name in seen:
            raise ValueError('unexpected or duplicate source path: ' + name)
        content = (item/'source'/name).read_bytes()
        if hashlib.sha256(content).hexdigest() != digest:
            raise ValueError('hash mismatch: ' + name)
        validate(name, content, pmcid)
        seen.add(name)
    if seen != set(assets) or {p.name for p in (item/'source').iterdir()} != seen:
        raise ValueError('source inventory differs from manifest')
    log = json.loads((item/'author_log.json').read_text())
    if log['dataset_id'] != dataset_id or {r['path'] for r in log['sources']} != seen:
        raise ValueError('author log source identity mismatch')
    for row in log['sources']:
        content = (item/'source'/row['path']).read_bytes()
        if row['sha256'] != hashlib.sha256(content).hexdigest() or row['bytes'] != len(content) or row['url'] != assets[row['path']]:
            raise ValueError('author log source provenance mismatch')
    print('PASS: ' + item.name + ' original assets, formats, hashes and provenance')


def acquire(item, dataset_id, pmcid, assets, provenance):
    """Download to temporary files and preserve pinned upstream bytes and metadata."""
    manifest = item / 'manifest.sha256'
    expected = {}
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            digest, name = line.split('  ', 1)
            expected[name] = digest
        if set(expected) != set(assets):
            raise ValueError('pinned manifest asset list differs')
    source = item / 'source'
    source.mkdir(exist_ok=True)
    records = []
    previous = json.loads((item/'author_log.json').read_text()) if (item/'author_log.json').exists() else None
    for name, url in assets.items():
        if Path(name).name != name:
            raise ValueError('source asset names must be basenames')
        dest = source / name
        if dest.exists():
            content = dest.read_bytes()
            media = next((r['media_type'] for r in previous['sources'] if r['path'] == name), 'unknown') if previous else 'unknown'
        else:
            with tempfile.TemporaryDirectory() as tmp:
                payload, headers = Path(tmp)/'payload', Path(tmp)/'headers'
                subprocess.run(['curl','-sS','-L','--fail','--max-time','30','--max-filesize','10000000','-D',str(headers),url,'-o',str(payload)], check=True)
                content = payload.read_bytes()
                media = next((line.split(':',1)[1].strip() for line in reversed(headers.read_text().splitlines()) if line.lower().startswith('content-type:')), 'unknown')
        if len(content) > 10_000_000:
            raise ValueError('unexpected asset size')
        validate(name, content, pmcid)
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
    if previous is None:
        log = dict(provenance, dataset_id=dataset_id, retrieved_at=datetime.now(timezone.utc).isoformat(), sources=records, local_additions=[])
        (item/'author_log.json').write_text(json.dumps(log,indent=2)+'\n')
    verify(item, dataset_id, pmcid, assets)
