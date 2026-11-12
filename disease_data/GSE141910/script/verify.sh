#!/usr/bin/env bash
# verify.sh — offline verification of GSE141910/source/ against ../manifest.sha256.
#
# The manifest lists `source/` only (ADR 0030), so `sha256sum -c` checks the
# downloaded bytes without touching local material and without the network.
# A missing manifest is a failure (exit 2), not a pass: regenerate it with
# `.scratch/data-provenance-restructure/ops/hash_items.py`.
set -uo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
MANIFEST="$ITEM_DIR/manifest.sha256"
if [ ! -f "$MANIFEST" ]; then
  echo "GSE141910: manifest missing — rerun .scratch/data-provenance-restructure/ops/hash_items.py" >&2
  exit 2
fi
cd -- "$ITEM_DIR/source" || { echo "GSE141910: source/ is missing" >&2; exit 2; }
echo "GSE141910: sha256sum -c ../manifest.sha256 (cwd: source/)"
sha256sum -c ../manifest.sha256
