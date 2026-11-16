#!/usr/bin/env bash
set -euo pipefail
ITEM="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ITEM/source"
if [ ! -s "$ITEM/source/GSE106297_RAW.tar" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106297/suppl/GSE106297_RAW.tar" -o "$ITEM/source/GSE106297_RAW.tar.part"
 mv "$ITEM/source/GSE106297_RAW.tar.part" "$ITEM/source/GSE106297_RAW.tar"
fi
if [ ! -s "$ITEM/source/GSE106297_series_matrix.txt.gz" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106297/matrix/GSE106297_series_matrix.txt.gz" -o "$ITEM/source/GSE106297_series_matrix.txt.gz.part"
 mv "$ITEM/source/GSE106297_series_matrix.txt.gz.part" "$ITEM/source/GSE106297_series_matrix.txt.gz"
fi
