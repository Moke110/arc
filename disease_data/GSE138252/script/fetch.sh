#!/usr/bin/env bash
# fetch.sh — download the GSE138252 deposit into source/ (ADR 0030).
#
# GSE138252 (PCM1+ FANS-sorted cardiomyocyte nuclei, bulk RNA-seq, Lipovsky et
# al. JCI Insight 2020) ships **no per-sample supplementary file**: every GSM's
# `!Sample_supplementary_file_1` is NONE and the data live only in the series'
# three supplementary matrices:
#   GSE138252_human_merged_genebody_exp.txt.gz                  1,114,989 B
#     56,609 genes x 18 human atrial libraries (the 2020 submission batch:
#     10 AF + 8 control); this is the item's primary payload.
#   GSE138252_merged_gene_name_expression.txt.gz                1,136,808 B
#     56,609 genes x 20 libraries (the 2021 addendum: 6 HF-only + 14 controls);
#     fetched and pinned, not used by the published contrast.
#   GSE138252_human_2785_merged_gene_name_expression.txt.gz       554,855 B
#     a 6-column subset of one sequencing run (technical replicates); fetched
#     and pinned, not used.
#   GSE138252_series_matrix.txt.gz                                  4,230 B
#     the !Sample_* design block (48 GSMs, tissue, disease state, molecule
#     subtype); its data table is empty.
# The deposits are tiny; a single-stream download is used.
#
# Usage: bash script/fetch.sh   (needs curl + network)
# Verify offline afterwards: bash script/verify.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
SERIES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138252"
mkdir -p "$SOURCE_DIR"

# url <TAB> upstream filename <TAB> expected bytes (GEO FTP listing)
FILES=(
  "$SERIES_BASE/matrix/GSE138252_series_matrix.txt.gz	GSE138252_series_matrix.txt.gz	4230"
  "$SERIES_BASE/suppl/GSE138252_human_merged_genebody_exp.txt.gz	GSE138252_human_merged_genebody_exp.txt.gz	1114989"
  "$SERIES_BASE/suppl/GSE138252_merged_gene_name_expression.txt.gz	GSE138252_merged_gene_name_expression.txt.gz	1136808"
  "$SERIES_BASE/suppl/GSE138252_human_2785_merged_gene_name_expression.txt.gz	GSE138252_human_2785_merged_gene_name_expression.txt.gz	554855"
)

download_one() {
  local url="$1" name="$2" bytes="$3" dest="$SOURCE_DIR/$2"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE138252: $name already present ($bytes B)"
    return 0
  fi
  echo "GSE138252: downloading $name"
  local attempt
  for attempt in 1 2 3; do
    if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --continue-at - --output "$dest.part" "$url"; then
      if [ "$(wc -c <"$dest.part" | tr -d ' ')" = "$bytes" ]; then
        mv -- "$dest.part" "$dest"
        echo "GSE138252: wrote source/$name ($bytes B)"
        return 0
      fi
    fi
    sleep $((attempt * 5))
  done
  echo "GSE138252: failed to fetch $name" >&2
  return 1
}

while IFS=$'\t' read -r url name bytes; do
  download_one "$url" "$name" "$bytes"
done < <(printf '%s\n' "${FILES[@]}")

while IFS=$'\t' read -r _url name bytes; do
  [ -f "$SOURCE_DIR/$name" ] || { echo "GSE138252: missing source/$name" >&2; exit 1; }
  [ "$(wc -c <"$SOURCE_DIR/$name" | tr -d ' ')" = "$bytes" ] \
    || { echo "GSE138252: wrong size for source/$name" >&2; exit 1; }
done < <(printf '%s\n' "${FILES[@]}")
echo "GSE138252: 4 payloads present; now run 'bash script/verify.sh'"
