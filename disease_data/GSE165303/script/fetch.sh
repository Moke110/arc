#!/usr/bin/env bash
# fetch.sh — download the GSE165303 deposit into source/ (ADR 0030).
#
# Two payloads only, both byte-size pinned below and re-checked after the
# download: the deposited RNA-seq master table (Ensembl gene x 101 libraries,
# log2 normalised counts plus the authors' DESeq2 columns) and the GPL20795
# series matrix (the deposit's own sample block, which binds those 101 columns to
# DCM/NF titles). A re-run re-downloads only a missing or wrong-sized file and
# never changes the bytes. The rest of the series (a 17.6 GB RAW.tar of HiChIP/
# Hi-C/ATAC payloads, and merged bigwigs) is deliberately not fetched: it carries
# no bulk RNA-seq count matrix. Derived expression/DEG products go to local/
# (declared in author_log.json).
#
# Usage: bash script/fetch.sh   (needs curl + network)
# Verify offline afterwards: bash script/verify.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
SERIES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165303"
mkdir -p "$SOURCE_DIR"

# url <TAB> upstream filename <TAB> expected bytes (GEO FTP listing)
FILES=(
  "$SERIES_BASE/suppl/GSE165303_RNA-Seq_master_table.txt.gz	GSE165303_RNA-Seq_master_table.txt.gz	20040757"
  "$SERIES_BASE/matrix/GSE165303-GPL20795_series_matrix.txt.gz	GSE165303-GPL20795_series_matrix.txt.gz	12562"
)

download_one() {
  local url="$1" name="$2" bytes="$3" dest="$SOURCE_DIR/$2"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE165303: $name already present ($bytes B)"
    return 0
  fi
  echo "GSE165303: downloading $name"
  local attempt
  for attempt in 1 2 3; do
    if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --continue-at - --output "$dest.part" "$url"; then
      if [ "$(wc -c <"$dest.part" | tr -d ' ')" = "$bytes" ]; then
        mv -- "$dest.part" "$dest"
        echo "GSE165303: wrote source/$name ($bytes B)"
        return 0
      fi
    fi
    sleep $((attempt * 5))
  done
  echo "GSE165303: failed to fetch $name" >&2
  return 1
}

while IFS=$'\t' read -r url name bytes; do
  download_one "$url" "$name" "$bytes"
done < <(printf '%s\n' "${FILES[@]}")

while IFS=$'\t' read -r _url name bytes; do
  [ -f "$SOURCE_DIR/$name" ] || { echo "GSE165303: missing source/$name" >&2; exit 1; }
  [ "$(wc -c <"$SOURCE_DIR/$name" | tr -d ' ')" = "$bytes" ] \
    || { echo "GSE165303: wrong size for source/$name" >&2; exit 1; }
done < <(printf '%s\n' "${FILES[@]}")
echo "GSE165303: 2 payloads present; now run 'bash script/verify.sh'"
