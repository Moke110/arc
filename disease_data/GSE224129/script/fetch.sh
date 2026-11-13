#!/usr/bin/env bash
# fetch.sh — download the GSE224129 deposit into source/ (ADR 0030).
#
# Two payloads only: the deposited expression matrix and the GEO series matrix
# (sample titles/characteristics — the genotype, cell-type, format and batch
# assignment of every library). Both are byte-size pinned below; a re-run
# re-downloads only a missing or wrong-sized file and never changes the bytes.
# Derived expression/DEG products go to local/ (declared in author_log.json).
#
# Usage: bash script/fetch.sh   (needs curl + network)
# Verify offline afterwards: bash script/verify.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
SERIES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224129"
mkdir -p "$SOURCE_DIR"

# upstream filename <TAB> expected bytes (GEO suppl/matrix listing)
FILES=(
  "GSE224129_AHNNF3DMXX_rsem_GENE.txt.gz 826308"
  "GSE224129_series_matrix.txt.gz 3834"
)

download_one() {
  local name="$1" bytes="$2" dest="$SOURCE_DIR/$1" url
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE224129: $name already present ($bytes B)"
    return 0
  fi
  case "$name" in
    *_series_matrix.txt.gz) url="$SERIES_BASE/matrix/$name" ;;
    *) url="$SERIES_BASE/suppl/$name" ;;
  esac
  echo "GSE224129: downloading $name"
  local attempt
  for attempt in 1 2 3; do
    if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --continue-at - --output "$dest.part" "$url"; then
      if [ "$(wc -c <"$dest.part" | tr -d ' ')" = "$bytes" ]; then
        mv -- "$dest.part" "$dest"
        echo "GSE224129: wrote source/$name ($bytes B)"
        return 0
      fi
    fi
    sleep $((attempt * 5))
  done
  echo "GSE224129: failed to fetch $name" >&2
  return 1
}

for spec in "${FILES[@]}"; do
  # shellcheck disable=SC2086
  download_one $spec
done
for spec in "${FILES[@]}"; do
  set -- $spec
  [ -f "$SOURCE_DIR/$1" ] || { echo "GSE224129: missing source/$1" >&2; exit 1; }
  [ "$(wc -c <"$SOURCE_DIR/$1" | tr -d ' ')" = "$2" ] \
    || { echo "GSE224129: wrong size for source/$1" >&2; exit 1; }
done
echo "GSE224129: 2 payloads present; now run 'bash script/verify.sh'"
