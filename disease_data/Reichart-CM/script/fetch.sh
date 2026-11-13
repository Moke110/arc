#!/usr/bin/env bash
# fetch.sh — download the cellxgene Reichart cardiomyocyte object into source/ (ADR 0030).
#
# This is a non-GEO item, named by its collection/object identity (as
# `pert_data/DToxS-iPSC` is). The object is one CM-only h5ad from the CELLxGENE
# collection "Pathogenic variants damage cell composition and single cell
# transcription in cardiomyopathies" (Reichart et al., Science 2022):
#
#   collection   e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb
#   dataset      "DCM/ACM heart cell atlas: Cardiomyocytes"
#   version      27fdbdf9-deae-4aed-bf7f-21964ee62137
#   asset        H5AD 3,194,863,189 B
#   url          https://datasets.cellxgene.cziscience.com/27fdbdf9-deae-4aed-bf7f-21964ee62137.h5ad
#   shape        311,418 cells; cell_type == "cardiac muscle cell" (CM-only object)
#
# The item keeps the downloaded h5ad under its version id so the object version
# is named in the filename as well as in author_log.json. source/ is the
# untouched download; the CM-gated per-donor pseudobulk, the genotype split and
# the tracking copies go to `local/` (declared in author_log.json).
#
# Usage: bash script/fetch.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
URL="https://datasets.cellxgene.cziscience.com/27fdbdf9-deae-4aed-bf7f-21964ee62137.h5ad"
NAME="27fdbdf9-deae-4aed-bf7f-21964ee62137.h5ad"
BYTES=3194863189
mkdir -p "$SOURCE_DIR"
dest="$SOURCE_DIR/$NAME"
if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$BYTES" ]; then
  echo "Reichart-CM: $NAME already present ($BYTES B)"
  exit 0
fi
echo "Reichart-CM: downloading $URL"
curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
  --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
  --continue-at - --output "$dest.part" "$URL"
actual=$(wc -c <"$dest.part" | tr -d ' ')
if [ "$actual" != "$BYTES" ]; then
  echo "Reichart-CM: $NAME is $actual B, expected $BYTES B" >&2
  exit 1
fi
mv -- "$dest.part" "$dest"
echo "Reichart-CM: wrote source/$NAME ($actual B)"
echo "Reichart-CM: now run 'bash script/verify.sh' (offline manifest check)"
