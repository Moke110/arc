#!/usr/bin/env bash
set -euo pipefail
ITEM="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ITEM/source"
if [ ! -s "$ITEM/source/GSE232331_hg38-counts.tsv.gz" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232331/suppl/GSE232331_hg38-counts.tsv.gz" -o "$ITEM/source/GSE232331_hg38-counts.tsv.gz.part"
 mv "$ITEM/source/GSE232331_hg38-counts.tsv.gz.part" "$ITEM/source/GSE232331_hg38-counts.tsv.gz"
fi
if [ ! -s "$ITEM/source/GSE232331_series_matrix.txt.gz" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232331/matrix/GSE232331_series_matrix.txt.gz" -o "$ITEM/source/GSE232331_series_matrix.txt.gz.part"
 mv "$ITEM/source/GSE232331_series_matrix.txt.gz.part" "$ITEM/source/GSE232331_series_matrix.txt.gz"
fi
