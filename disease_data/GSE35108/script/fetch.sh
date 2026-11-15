#!/usr/bin/env bash
# fetch.sh — download the GSE35108 deposit into source/ (ADR 0030).
#
# Two payloads only, both byte-size pinned below and re-checked after the
# download: the GEO series matrix (the deposited RMA log2 expression table plus
# the !Sample_* design block) and the GPL6244 platform annotation (the
# Affymetrix transcript-cluster -> Entrez gene symbol table this one-colour array
# needs). A re-run re-downloads only a missing or wrong-sized file and never
# changes the bytes. Derived expression/DEG products go to local/ (declared in
# author_log.json).
#
# Usage: bash script/fetch.sh   (needs curl + network)
# Verify offline afterwards: bash script/verify.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
SERIES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE35nnn/GSE35108"
PLATFORM_BASE="https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244"
mkdir -p "$SOURCE_DIR"

# url <TAB> upstream filename <TAB> expected bytes (GEO FTP listing)
FILES=(
  "$SERIES_BASE/matrix/GSE35108_series_matrix.txt.gz	GSE35108_series_matrix.txt.gz	2302867"
  "$PLATFORM_BASE/annot/GPL6244.annot.gz	GPL6244.annot.gz	7193437"
)

download_one() {
  local url="$1" name="$2" bytes="$3" dest="$SOURCE_DIR/$2"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE35108: $name already present ($bytes B)"
    return 0
  fi
  echo "GSE35108: downloading $name"
  local attempt
  for attempt in 1 2 3; do
    if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --continue-at - --output "$dest.part" "$url"; then
      if [ "$(wc -c <"$dest.part" | tr -d ' ')" = "$bytes" ]; then
        mv -- "$dest.part" "$dest"
        echo "GSE35108: wrote source/$name ($bytes B)"
        return 0
      fi
    fi
    sleep $((attempt * 5))
  done
  echo "GSE35108: failed to fetch $name" >&2
  return 1
}

while IFS=$'\t' read -r url name bytes; do
  download_one "$url" "$name" "$bytes"
done < <(printf '%s\n' "${FILES[@]}")

while IFS=$'\t' read -r _url name bytes; do
  [ -f "$SOURCE_DIR/$name" ] || { echo "GSE35108: missing source/$name" >&2; exit 1; }
  [ "$(wc -c <"$SOURCE_DIR/$name" | tr -d ' ')" = "$bytes" ] \
    || { echo "GSE35108: wrong size for source/$name" >&2; exit 1; }
done < <(printf '%s\n' "${FILES[@]}")
echo "GSE35108: 2 payloads present; now run 'bash script/verify.sh'"
