#!/usr/bin/env bash
# Verify the downloaded HGNC source (ADR 0030).
# A missing manifest is a failure (exit 2), not a pass: regenerate it with
# `.scratch/data-provenance-restructure/ops/hash_items.py`.
set -uo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
MANIFEST="$ITEM_DIR/manifest.sha256"
if [ ! -f "$MANIFEST" ]; then
  echo "Annotation: manifest missing — rerun .scratch/data-provenance-restructure/ops/hash_items.py" >&2
  exit 2
fi
cd -- "$ITEM_DIR/source" || { echo "Annotation: source/ is missing" >&2; exit 2; }
echo "Annotation: sha256sum -c ../manifest.sha256 (cwd: source/)"
sha256sum -c ../manifest.sha256
