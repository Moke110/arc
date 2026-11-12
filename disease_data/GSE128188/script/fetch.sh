#!/usr/bin/env bash
# fetch.sh — download the GSE128188 upstream supplement into source/ (ADR 0030).
#
# The item keeps a single upstream file, byte-identical to GEO's copy:
#
#   GSE128188_CountsV2_EdgeR.txt.gz   647,785 B
#   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128188/suppl/
#
# `source/` is the untouched download (ADR 0030); the derived DEG/expression
# tables are written by script/build_deg.R and script/build_expression.R into
# `local/`, never here. Re-running this script replaces source/ with a fresh
# upstream copy; compare against manifest.sha256 to prove byte identity.
#
# Usage:
#   bash script/fetch.sh
#
# No R, no package installs; needs curl (or wget) and the network. The offline
# check afterwards is `bash script/verify.sh`.

set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128188/suppl"
FILE="GSE128188_CountsV2_EdgeR.txt.gz"
EXPECTED_BYTES=647785

mkdir -p "$SOURCE_DIR"
dest="$SOURCE_DIR/$FILE"

echo "GSE128188: downloading $BASE_URL/$FILE"
if command -v curl >/dev/null 2>&1; then
  curl --fail --silent --show-error --location --retry 3 \
    --output "$dest.part" "$BASE_URL/$FILE"
else
  wget --quiet --tries=3 --output-document="$dest.part" "$BASE_URL/$FILE"
fi
mv -- "$dest.part" "$dest"

actual_bytes=$(wc -c <"$dest" | tr -d ' ')
if [ "$actual_bytes" -ne "$EXPECTED_BYTES" ]; then
  echo "GSE128188: $FILE is $actual_bytes B, expected $EXPECTED_BYTES B" >&2
  exit 1
fi
echo "GSE128188: wrote source/$FILE ($actual_bytes B)"
echo "GSE128188: now run 'bash script/verify.sh' (offline manifest check)"
