#!/usr/bin/env bash
# fetch.sh — download the GSE255612 left-atrium snRNA deposit into source/ (ADR 0030).
#
# GSE255612 (human left atrium, 10x 3' v3 snRNA) ships its usable payload as one
# Scanpy Matrix-Market triplet plus the author cell table:
#   GSE255612_AF_snRNA_Matrix_V1.mtx.gz                                     1,003,414,933 B
#     feature x cell, 36,600 x 179,697, 288,701,961 nonzeros; REAL values. The
#     series-matrix processing note says "Processed data files output from
#     Scanpy as gene expression matrix (.mtx)"; the values are log1p(CP10K)
#     (a cell's single-count genes share one value, = log1p(1e4/n_umi)).
#   GSE255612_AF_snRNA_Processed_Expression_Matrix_genes_final.tsv.gz        277,568 B
#     the mtx feature axis (Ensembl id + symbol); without it the matrix rows are
#     unnameable, so it is part of the deposit, not an extra payload.
#   GSE255612_AF_snRNA_MetaData.txt.gz                                     8,846,174 B
#     the cell table (NAME barcode, donor_id, cell_type, af case/control, ...);
#     its 179,697 rows are in mtx column order (verified against the deposited
#     barcodes_V1 file, which is NOT stored: barcode order is reproducible from
#     the metadata itself).
# The series RAW.tar is empty for this deposit; no other file is fetched.
#
# source/ is the untouched download; the CM-gated pseudobulk goes to `local/`
# (declared in author_log.json). GEO throttles a single connection, so the two
# large files are fetched as parallel byte ranges and concatenated; the bytes are
# identical to a single-stream download and are checked by size + gzip integrity.
#
# Usage: bash script/fetch.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255612/suppl"
WORKERS=${GSE255612_WORKERS:-8}
mkdir -p "$SOURCE_DIR"

# Parallel ranged download of one URL; curl per range, then cat in order.
download_ranges() {
  local name="$1" bytes="$2"
  local dest="$SOURCE_DIR/$name"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE255612: $name already present ($bytes B)"
    return 0
  fi
  local tmp="$SOURCE_DIR/.$name.parts"
  rm -rf -- "$tmp"; mkdir -p "$tmp"
  local chunk=$(( (bytes + WORKERS - 1) / WORKERS ))
  local i=0
  while [ $(( i * chunk )) -lt "$bytes" ]; do
    local start=$(( i * chunk ))
    local end=$(( start + chunk - 1 ))
    [ "$end" -ge "$bytes" ] && end=$(( bytes - 1 ))
    echo "GSE255612: $name range $start-$end"
    (curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --range "$start-$end" --output "$tmp/part.$i" "$BASE/$name" \
      || { echo "GSE255612: range $i failed" >&2; exit 1; }) &
    i=$(( i + 1 ))
    if [ $(( i % WORKERS )) -eq 0 ]; then wait; fi
  done
  wait
  # Concatenate in range order, never through a shell glob: `part.10` sorts
  # before `part.2`, so a lexicographic glob corrupts any run with more than
  # nine ranges (the size check still passes; only gzip -t catches it). The
  # explicit index loop also fails closed when a range never landed, and a
  # failed run leaves no part file behind for source/ to hash.
  local n_parts="$i"
  local part=0
  local part_file
  : >"$dest.part"
  while [ "$part" -lt "$n_parts" ]; do
    part_file="$tmp/part.$part"
    if [ ! -f "$part_file" ]; then
      echo "GSE255612: missing range part $part" >&2
      rm -rf -- "$tmp" "$dest.part"
      exit 1
    fi
    cat -- "$part_file" >>"$dest.part"
    part=$(( part + 1 ))
  done
  rm -rf -- "$tmp"
  local actual
  actual=$(wc -c <"$dest.part" | tr -d ' ')
  if [ "$actual" != "$bytes" ]; then
    echo "GSE255612: $name is $actual B, expected $bytes B" >&2
    rm -f -- "$dest.part"
    exit 1
  fi
  mv -- "$dest.part" "$dest"
  echo "GSE255612: wrote source/$name ($bytes B)"
}

download_ranges "GSE255612_AF_snRNA_Matrix_V1.mtx.gz" 1003414933
download_ranges "GSE255612_AF_snRNA_Processed_Expression_Matrix_genes_final.tsv.gz" 277568

meta="GSE255612_AF_snRNA_MetaData.txt.gz"
if [ ! -f "$SOURCE_DIR/$meta" ] || [ "$(wc -c <"$SOURCE_DIR/$meta" | tr -d ' ')" != 8846174 ]; then
  echo "GSE255612: downloading $meta"
  curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
    --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
    --continue-at - --output "$SOURCE_DIR/$meta.part" "$BASE/$meta"
  mv -- "$SOURCE_DIR/$meta.part" "$SOURCE_DIR/$meta"
fi
echo "GSE255612: wrote source/$meta"

# gzip integrity is a cheap byte-level guard on the ranged download.
for f in "$SOURCE_DIR"/*.gz; do
  gzip -t -- "$f" || { echo "GSE255612: corrupted gzip: $f" >&2; exit 1; }
done
echo "GSE255612: now run 'bash script/verify.sh' (offline manifest check)"
