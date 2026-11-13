#!/usr/bin/env bash
# verify.sh — offline verification of GSE270788/source/ against ../manifest.sha256.
#
# The manifest lists `source/` only (ADR 0030), so `sha256sum -c` checks the
# downloaded bytes without touching local material and without the network.
# A missing manifest is a failure (exit 2), not a pass: regenerate it with
# script/hash_source.py (tracked, ADR 0030 rule: sorted source-relative paths).
set -uo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
MANIFEST="$ITEM_DIR/manifest.sha256"
if [ ! -s "$MANIFEST" ]; then
  echo "GSE270788: manifest.sha256 missing — regenerate with script/hash_source.py" >&2
  exit 2
fi
cd -- "$ITEM_DIR/source" || { echo "GSE270788: source/ is missing" >&2; exit 2; }
echo "GSE270788: sha256sum -c ../manifest.sha256 (cwd: source/)"
sha256sum -c ../manifest.sha256
