#!/usr/bin/env python
"""hash_source.py — regenerate an ADR 0030 `manifest.sha256` over `source/`.

The ADR 0030 manifest rule (recorded at acceptance, 2026-10-30) is one
`sha256sum` line per file under `source/`, sorted by **source-relative path**;
`tree_sha256` in `author_log.json` is sha256 over the manifest file bytes (the
clarification explicitly rejects a byte-wise lexicographic sort of the rendered
lines). The historical generator `.scratch/.../hash_items.py` no longer exists,
so this tracked copy is the rule's live implementation. It is run as:

    /home/chang/miniconda3/envs/arc/bin/python script/hash_source.py

The result is written to `<item>/manifest.sha256`; verification is the offline
`cd source && sha256sum -c ../manifest.sha256` in script/verify.sh.
"""
import hashlib
import sys
from pathlib import Path

ITEM = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    source = ITEM / "source"
    files = sorted((p for p in source.rglob("*") if p.is_file()),
                   key=lambda p: p.relative_to(source).as_posix())
    if not files:
        raise SystemExit(f"{ITEM.name}: source/ is empty")
    lines = [f"{sha256(path)}  {path.relative_to(source).as_posix()}" for path in files]
    manifest = ITEM / "manifest.sha256"
    payload = ("\n".join(lines) + "\n").encode()
    manifest.write_bytes(payload)
    print(f"{ITEM.name}: {len(lines)} file(s), manifest.sha256")
    print(f"{ITEM.name}: tree_sha256 {hashlib.sha256(payload).hexdigest()}")
    print(f"{ITEM.name}: source bytes {sum(p.stat().st_size for p in files)}")


if __name__ == "__main__":
    sys.exit(main())
