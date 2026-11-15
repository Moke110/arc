#!/usr/bin/env bash
# fetch.sh — download the GSE238242 snRNA arm into source/ (ADR 0030).
#
# GSE238242 is a 10x multiome deposit (paired snATAC + snRNA, 14 libraries) plus
# an author cell-type table. The series' usable deposit is `GSE238242_RAW.tar`
# (187,269,120 B = 7 RNA + 7 ATAC dense count tables); ADR 0071's source is the
# **snRNA arm only**, so the seven RNA libraries are fetched per-GSM from GEO's
# per-sample supplementary directory (~45 MB) and the ~130 MB ATAC arm is never
# downloaded. Every fetched byte is checked against the size in the series'
# GEO file listing (captured in `author_log.json` notes).
#
#   GSE238242_series_matrix.txt.gz                               3,276 B
#     the !Sample_* design block (14 libraries, titles with rhythm and donor,
#     `rhythm: SR/AF`); its data table is empty.
#   GSE238242_snAF.metadata.tsv.gz                              74,550 B
#     the author cell table: barcode -> cell_type / sample / sex / Rhythm;
#     `cell_type == CM` gates the cardiomyocyte nuclei.
#   GSM7660990_CF69RNA.counts.tsv.gz ... GSM7660996_CF102RNA.counts.tsv.gz
#     the seven RNA count tables (dense TSV, genes x nuclei, raw integers;
#     38,231 genes each), byte-identical to the tar members.
#
# Usage: bash script/fetch.sh   (needs curl + network)
# Verify offline afterwards: bash script/verify.sh
set -euo pipefail
ITEM_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$ITEM_DIR/source"
SERIES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238242"
SAMPLE_BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7660nnn"
mkdir -p "$SOURCE_DIR"

# url <TAB> upstream filename <TAB> expected bytes
FILES=(
  "$SERIES_BASE/matrix/GSE238242_series_matrix.txt.gz	GSE238242_series_matrix.txt.gz	3276"
  "$SERIES_BASE/suppl/GSE238242_snAF.metadata.tsv.gz	GSE238242_snAF.metadata.tsv.gz	74550"
  "$SAMPLE_BASE/GSM7660990/suppl/GSM7660990_CF69RNA.counts.tsv.gz	GSM7660990_CF69RNA.counts.tsv.gz	9441961"
  "$SAMPLE_BASE/GSM7660991/suppl/GSM7660991_CF77RNA.counts.tsv.gz	GSM7660991_CF77RNA.counts.tsv.gz	7136770"
  "$SAMPLE_BASE/GSM7660992/suppl/GSM7660992_CF89RNA.counts.tsv.gz	GSM7660992_CF89RNA.counts.tsv.gz	3547052"
  "$SAMPLE_BASE/GSM7660993/suppl/GSM7660993_CF91RNA.counts.tsv.gz	GSM7660993_CF91RNA.counts.tsv.gz	5856219"
  "$SAMPLE_BASE/GSM7660994/suppl/GSM7660994_CF93RNA.counts.tsv.gz	GSM7660994_CF93RNA.counts.tsv.gz	5468856"
  "$SAMPLE_BASE/GSM7660995/suppl/GSM7660995_CF97RNA.counts.tsv.gz	GSM7660995_CF97RNA.counts.tsv.gz	8797614"
  "$SAMPLE_BASE/GSM7660996/suppl/GSM7660996_CF102RNA.counts.tsv.gz	GSM7660996_CF102RNA.counts.tsv.gz	5107455"
)

download_one() {
  local url="$1" name="$2" bytes="$3" dest="$SOURCE_DIR/$2"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest" | tr -d ' ')" = "$bytes" ]; then
    echo "GSE238242: $name already present ($bytes B)"
    return 0
  fi
  echo "GSE238242: downloading $name"
  local attempt
  for attempt in 1 2 3; do
    if curl --fail --silent --show-error --location --retry 5 --retry-delay 5 \
        --connect-timeout 30 --speed-limit 2048 --speed-time 60 \
        --continue-at - --output "$dest.part" "$url"; then
      if [ "$(wc -c <"$dest.part" | tr -d ' ')" = "$bytes" ]; then
        mv -- "$dest.part" "$dest"
        echo "GSE238242: wrote source/$name ($bytes B)"
        return 0
      fi
    fi
    sleep $((attempt * 5))
  done
  echo "GSE238242: failed to fetch $name" >&2
  return 1
}

while IFS=$'\t' read -r url name bytes; do
  download_one "$url" "$name" "$bytes"
done < <(printf '%s\n' "${FILES[@]}")

for f in "$SOURCE_DIR"/*.gz; do
  gzip -t -- "$f" || { echo "GSE238242: corrupted gzip: $f" >&2; exit 1; }
done
while IFS=$'\t' read -r _url name bytes; do
  [ -f "$SOURCE_DIR/$name" ] || { echo "GSE238242: missing source/$name" >&2; exit 1; }
  [ "$(wc -c <"$SOURCE_DIR/$name" | tr -d ' ')" = "$bytes" ] \
    || { echo "GSE238242: wrong size for source/$name" >&2; exit 1; }
done < <(printf '%s\n' "${FILES[@]}")
echo "GSE238242: 9 payloads present (snRNA arm only; ATAC arm not fetched); now run 'bash script/verify.sh'"
