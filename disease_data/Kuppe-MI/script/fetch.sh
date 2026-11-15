#!/usr/bin/env bash
# fetch.sh — download the cellxgene Kuppe MI all-snRNA object into source/ (ADR 0030).
#
# This is a non-GEO item, named by its collection/object identity (as
# `pert_data/DToxS-iPSC` and `disease_data/Reichart-CM` are). The object is one
# h5ad from the CELLxGENE collection "Spatial multi-omic map of human myocardial
# infarction" (Kuppe et al., Nature 2022):
#
#   collection   8191c283-0816-424b-9b61-c3e1d6258a77
#   dataset      "All-snRNA-Spatial multi-omic map of human myocardial infarction"
#                dataset id 1c739a3e-c3f5-49d5-98e0-73975e751201
#   version      4780eb7d-4c7c-4630-8705-471f0de8e27a
#   asset        H5AD 1,862,945,764 B
#   url          https://datasets.cellxgene.cziscience.com/4780eb7d-4c7c-4630-8705-471f0de8e27a.h5ad
#   shape        191,795 nuclei, 20 donors, 28,380 features, left ventricle
#                cell_type_original == "Cardiomyocyte" (64,510 nuclei)
#
# The version id is kept in the filename so the object version is named there as
# well as in author_log.json. `source/` is the untouched download; the CM-gated
# per-(donor, zone) pseudobulk and the DEG products go to `local/` (declared in
# author_log.json's local_additions[]). The canonical recipe is exactly this
# curl: the cellxgene asset endpoint serves the version's bytes; a re-run
# re-downloads only a missing or wrong-sized file and never changes the bytes.
#
# Usage: bash script/fetch.sh   (needs curl + network)
# Verify offline afterwards: bash script/verify.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
URL="https://datasets.cellxgene.cziscience.com/4780eb7d-4c7c-4630-8705-471f0de8e27a.h5ad"
NAME="4780eb7d-4c7c-4630-8705-471f0de8e27a.h5ad"
BYTES=1862945764
mkdir -p "$SOURCE_DIR"
dest="$SOURCE_DIR/$NAME"
if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$BYTES" ]; then
  echo "Kuppe-MI: $NAME already present ($BYTES B)"
  exit 0
fi
echo "Kuppe-MI: downloading $URL"
attempt=0
for attempt in 1 2 3; do
  if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
      --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
      --continue-at - --output "$dest.part" "$URL"; then
    actual=$(wc -c <"$dest.part" | tr -d ' ')
    if [ "$actual" = "$BYTES" ]; then
      mv -- "$dest.part" "$dest"
      echo "Kuppe-MI: wrote source/$NAME ($actual B)"
      echo "Kuppe-MI: now run 'bash script/verify.sh' (offline manifest check)"
      exit 0
    fi
    echo "Kuppe-MI: $NAME is $actual B, expected $BYTES B" >&2
  fi
  sleep $((attempt * 5))
done
echo "Kuppe-MI: failed to fetch $NAME" >&2
exit 1
