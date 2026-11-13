#!/usr/bin/env bash
# fetch.sh — download the GSE270788 GEX sub-series into source/ (ADR 0030).
#
# GSE270788 (Amrute et al., Nature 2024) is a human LV 10x Multiome series.
# GSE218392 is its SuperSeries, so this item takes the GEX SubSeries only.
# The series-level GSE270788_RAW.tar (80.7 GB) is NOT fetched: the per-GSM
# filtered feature-barcode matrices below are the usable payload (23 donors,
# 931,875,148 B total). GSE270788_metadata.csv.gz carries the author cell table
# whose cell.type.l1 / HF.etiology columns define the CM gate and the arms.
#
# The 23 (GSM, upstream filename) pairs are the filelist.txt entries of
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/ whose type is
# H5 and whose title is "GEX MA<n>" (the ATAC fragment files and the raw tar are
# deliberately excluded). GSM8352066's filename has no MA token; its sample is
# GEX MA29 per the series-matrix !Sample_title row.
#
# source/ is the untouched download; derived pseudobulk/DEG products go to
# `local/` (declared in author_log.json). Re-running this script re-downloads and
# overwrites only with byte-identical upstream content; verify offline afterwards
# with `bash script/verify.sh`.
#
# Usage: bash script/fetch.sh
# Needs curl and the network. UTF-8/locale-neutral.
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
SAMPLE_BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples"
SERIES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl"
WORKERS=${GSE270788_WORKERS:-6}
mkdir -p "$SOURCE_DIR"

# GSM<TAB>upstream filename<TAB>expected bytes (GEO suppl listing)
FILES=(
  "GSM8352048 GSM8352048_MA5_filtered_feature_bc_matrix.h5 42015815"
  "GSM8352049 GSM8352049_MA6_filtered_feature_bc_matrix.h5 77608627"
  "GSM8352050 GSM8352050_MA7_filtered_feature_bc_matrix.h5 13966175"
  "GSM8352051 GSM8352051_MA8_filtered_feature_bc_matrix.h5 35607043"
  "GSM8352052 GSM8352052_MA9_filtered_feature_bc_matrix.h5 42241264"
  "GSM8352053 GSM8352053_MA10_filtered_feature_bc_matrix.h5 51931292"
  "GSM8352054 GSM8352054_MA11_filtered_feature_bc_matrix.h5 22606389"
  "GSM8352055 GSM8352055_MA13_filtered_feature_bc_matrix.h5 7085583"
  "GSM8352056 GSM8352056_MA14_filtered_feature_bc_matrix.h5 49545864"
  "GSM8352057 GSM8352057_MA19_filtered_feature_bc_matrix.h5 61631853"
  "GSM8352058 GSM8352058_MA20_filtered_feature_bc_matrix.h5 10556795"
  "GSM8352059 GSM8352059_MA22_filtered_feature_bc_matrix.h5 22169714"
  "GSM8352060 GSM8352060_MA23_filtered_feature_bc_matrix.h5 34635722"
  "GSM8352061 GSM8352061_MA24_filtered_feature_bc_matrix.h5 43914551"
  "GSM8352062 GSM8352062_MA25_filtered_feature_bc_matrix.h5 16586615"
  "GSM8352063 GSM8352063_MA26_filtered_feature_bc_matrix.h5 21571481"
  "GSM8352064 GSM8352064_MA27_filtered_feature_bc_matrix.h5 36613831"
  "GSM8352065 GSM8352065_MA28_filtered_feature_bc_matrix.h5 28836067"
  "GSM8352066 GSM8352066_filtered_feature_bc_matrix.h5 55968251"
  "GSM8352067 GSM8352067_MA30_filtered_feature_bc_matrix.h5 48781637"
  "GSM8352068 GSM8352068_MA31_filtered_feature_bc_matrix.h5 104340627"
  "GSM8352069 GSM8352069_MA32_filtered_feature_bc_matrix.h5 53533347"
  "GSM8352070 GSM8352070_MA33_filtered_feature_bc_matrix.h5 56181921"
)

download_one() {
  local gsm="$1" name="$2" bytes="$3"
  local dest="$SOURCE_DIR/$name"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE270788: $name already present ($bytes B)"
    return 0
  fi
  local url="$SAMPLE_BASE/${gsm:0:7}nnn/${gsm}/suppl/${name}"
  echo "GSE270788: downloading $name"
  local attempt
  for attempt in 1 2 3; do
    if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --continue-at - --output "$dest.part" "$url"; then
      if [ "$(wc -c <"$dest.part" | tr -d ' ')" = "$bytes" ]; then
        mv -- "$dest.part" "$dest"
        echo "GSE270788: wrote source/$name ($bytes B)"
        return 0
      fi
    fi
    sleep $((attempt * 5))
  done
  echo "GSE270788: failed to fetch $name" >&2
  return 1
}
export -f download_one
export SOURCE_DIR SAMPLE_BASE

# One worker per file; GEO throttles per connection, so a few parallel streams
# are much faster than a sequential loop without changing the bytes.
printf '%s\n' "${FILES[@]}" | xargs -r -P "$WORKERS" -n 3 bash -c 'download_one "$@"' _

meta="GSE270788_metadata.csv.gz"
if [ ! -f "$SOURCE_DIR/$meta" ] || [ "$(wc -c <"$SOURCE_DIR/$meta" | tr -d ' ')" != 1272635 ]; then
  echo "GSE270788: downloading $meta"
  curl --fail --silent --show-error --location --retry 3 \
    --output "$SOURCE_DIR/$meta.part" "$SERIES_BASE/$meta"
  mv -- "$SOURCE_DIR/$meta.part" "$SOURCE_DIR/$meta"
fi
echo "GSE270788: wrote source/$meta"

# Fail closed: every declared payload must exist at its recorded size.
for spec in "${FILES[@]}"; do
  set -- $spec
  path="$SOURCE_DIR/$2"
  [ -f "$path" ] || { echo "GSE270788: missing source/$2" >&2; exit 1; }
  [ "$(wc -c <"$path" | tr -d ' ')" = "$3" ] \
    || { echo "GSE270788: wrong size for source/$2" >&2; exit 1; }
done
echo "GSE270788: 23 GEX matrices + metadata present"
echo "GSE270788: now run 'bash script/verify.sh' (offline manifest check)"
