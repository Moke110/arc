#!/usr/bin/env bash
set -euo pipefail
ITEM="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ITEM/source"
if [ ! -s "$ITEM/source/GSE235470_rawcounts.counts.txt.gz" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235470/suppl/GSE235470_rawcounts.counts.txt.gz" -o "$ITEM/source/GSE235470_rawcounts.counts.txt.gz.part"
 mv "$ITEM/source/GSE235470_rawcounts.counts.txt.gz.part" "$ITEM/source/GSE235470_rawcounts.counts.txt.gz"
fi
if [ ! -s "$ITEM/source/GSE235470_series_matrix.txt.gz" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235470/matrix/GSE235470_series_matrix.txt.gz" -o "$ITEM/source/GSE235470_series_matrix.txt.gz.part"
 mv "$ITEM/source/GSE235470_series_matrix.txt.gz.part" "$ITEM/source/GSE235470_series_matrix.txt.gz"
fi
