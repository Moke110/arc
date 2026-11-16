#!/usr/bin/env bash
set -euo pipefail
ITEM="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ITEM/source"
fetch() {
 local url="$1" name="$2"
 if [ ! -s "$ITEM/source/$name" ]; then
  curl -fL --retry 3 "$url" -o "$ITEM/source/$name.part"
  mv "$ITEM/source/$name.part" "$ITEM/source/$name"
 fi
}
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313746/suppl/GSE313746_Counts_RNA_EMPfortmiller.txt.gz" "GSE313746_Counts_RNA_EMPfortmiller.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313746/matrix/GSE313746_series_matrix.txt.gz" "GSE313746_series_matrix.txt.gz"
