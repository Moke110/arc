#!/usr/bin/env bash
# verify.sh — offline verification of GSE244740/source/ against ../manifest.sha256.
#
# The manifest lists `source/` only (ADR 0030), so `sha256sum -c` checks the
# downloaded bytes without touching local material and without the network.
# The manifest is produced by ticket 09; until then this is a declared no-op.
set -uo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
MANIFEST="$ITEM_DIR/manifest.sha256"
if [ ! -f "$MANIFEST" ]; then
  echo "GSE244740: manifest pending (ticket 09) — nothing to verify yet"
  exit 0
fi
cd -- "$ITEM_DIR/source" || { echo "GSE244740: source/ is missing" >&2; exit 2; }
echo "GSE244740: sha256sum -c ../manifest.sha256 (cwd: source/)"
sha256sum -c ../manifest.sha256
