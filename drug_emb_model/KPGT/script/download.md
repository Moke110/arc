# Download — KPGT

**method: vendored copy (manual)** — no upstream URL or commit is reproducible
from this host, so `author_log.json` records `version: "unknown"` and the four
later-mtime files below are of unknown divergence.

## What is here

- `source/{src,scripts}/`, `source/environment.yml`, `source/LICENSE`,
  `source/.python-version`, `source/.gitignore`, `source/KPGT.png` — the KPGT
  code as vendored (Apache-2.0, see `source/LICENSE`).
- `source/models/pretrained/base/{base.pth,base_cpu.pth}` — downloaded pretrained
  base weights, 894,215,942 B total; declared as an `asset` in `author_log.json`.
- `source/feature_extract/{LINCS2020,GSE92742,FDA_drugs}/` — downloaded input
  feature tensors/tables (34,997,345 B); declared as an `asset`.

## Divergence note

All code files are mtime 2026-04-18 14:46 except four whose mtime is later and
whose divergence from upstream is **unknown** (there is no upstream reference
tree on this host):

- `source/src/utils.py` — 2026-04-25 00:27
- `source/scripts/extract_features.py` — 2026-04-25 00:28
- `source/src/smiles_to_kpgt_pickle.py` — 2026-04-25 00:28
- `source/src/canonicalize_smiles.py` — 2026-04-25 00:44

Treat those four as possibly locally modified; re-derive from upstream before
trusting byte-level reproducibility.
