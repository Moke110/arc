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
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242692/suppl/GSE242692_countsMatrix_raw.txt.gz" "GSE242692_countsMatrix_raw.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242692/matrix/GSE242692_series_matrix.txt.gz" "GSE242692_series_matrix.txt.gz"
