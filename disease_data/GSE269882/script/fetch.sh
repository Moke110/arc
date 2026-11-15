#!/usr/bin/env bash
# fetch.sh — download the GSE269882 deposit into source/ (ADR 0030).
#
# GSE269882 (right atrial appendage bulk RNA-seq; SR / post-operative AF /
# paroxysmal AF / persistent AF, 27 libraries) ships three series-level files
# and no per-sample supplementary data:
#   GSE269882_countdata_all.csv.gz                             1,352,166 B
#     raw featureCounts integers, 57,500 genes x 27 samples (genes as
#     Ensembl ids in column 1, sample columns named by patient number);
#     the item's DEG payload.
#   GSE269882_TPM_values.csv.gz                                2,875,477 B
#     the deposited TPM matrix for the same 27 samples; fetched and pinned
#     for the cross-check, not used as the DEG input (the count model is
#     the documented one).
#   GSE269882_series_matrix.txt.gz                                 3,636 B
#     the !Sample_* design block: titles with patient numbers, Sex, age and
#     `diagnosis group`; its data table is empty.
#
# Usage: bash script/fetch.sh   (needs curl + network)
# Verify offline afterwards: bash script/verify.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
SERIES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269882"
mkdir -p "$SOURCE_DIR"

FILES=(
  "$SERIES_BASE/matrix/GSE269882_series_matrix.txt.gz	GSE269882_series_matrix.txt.gz	3636"
  "$SERIES_BASE/suppl/GSE269882_countdata_all.csv.gz	GSE269882_countdata_all.csv.gz	1352166"
  "$SERIES_BASE/suppl/GSE269882_TPM_values.csv.gz	GSE269882_TPM_values.csv.gz	2875477"
)

download_one() {
  local url="$1" name="$2" bytes="$3" dest="$SOURCE_DIR/$2"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE269882: $name already present ($bytes B)"
    return 0
  fi
  echo "GSE269882: downloading $name"
  local attempt
  for attempt in 1 2 3; do
    if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --continue-at - --output "$dest.part" "$url"; then
      if [ "$(wc -c <"$dest.part" | tr -d ' ')" = "$bytes" ]; then
        mv -- "$dest.part" "$dest"
        echo "GSE269882: wrote source/$name ($bytes B)"
        return 0
      fi
    fi
    sleep $((attempt * 5))
  done
  echo "GSE269882: failed to fetch $name" >&2
  return 1
}

while IFS=$'\t' read -r url name bytes; do
  download_one "$url" "$name" "$bytes"
done < <(printf '%s\n' "${FILES[@]}")

while IFS=$'\t' read -r _url name bytes; do
  [ -f "$SOURCE_DIR/$name" ] || { echo "GSE269882: missing source/$name" >&2; exit 1; }
  [ "$(wc -c <"$SOURCE_DIR/$name" | tr -d ' ')" = "$bytes" ] \
    || { echo "GSE269882: wrong size for source/$name" >&2; exit 1; }
done < <(printf '%s\n' "${FILES[@]}")
echo "GSE269882: 3 payloads present; now run 'bash script/verify.sh'"
