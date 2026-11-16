#!/usr/bin/env bash
# download.sh — fetch the GSE244740 GEO submission into source/ (ADR 0030).
#
# method: dataset — NCBI GEO public FTP; no authenticated or browser step.
# source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244740
# The two file names are the ones actually present in source/.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ITEM_DIR="$(dirname -- "$HERE")"
ROOT="$ITEM_DIR"
while [ ! -e "$ROOT/.git" ]; do
  [ "$ROOT" = "/" ] && { echo "download: cannot find the repository root" >&2; exit 2; }
  ROOT="$(dirname -- "$ROOT")"
done
GSE=GSE244740
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/${GSE}"
OUT="$ITEM_DIR/source"
mkdir -p "$OUT"

fetch() { # fetch <url> <name>
  local url="$1" name="$2"
  if [ -s "$OUT/$name" ]; then echo "have  $name"; return 0; fi
  echo "fetch $name"
  curl -fL --retry 3 --retry-delay 2 -o "$OUT/$name.part" "$url"
  mv -- "$OUT/$name.part" "$OUT/$name"
}

# processed counts: GEO supplementary file (.../suppl/)
fetch "${BASE}/suppl/${GSE}_processed_data_counts.txt.gz" "${GSE}_processed_data_counts.txt.gz"
# series matrix: GEO matrix file (.../matrix/)
fetch "${BASE}/matrix/${GSE}_series_matrix.txt.gz" "${GSE}_series_matrix.txt.gz"
# Author-hosted primary paper: Table 1 independently labels the reference panel.
fetch "https://news.cyprotex.com/hubfs/Downloads/Articles%20and%20white%20papers/In-depth%20mechanistic%20analysis%20including%20high-throughput%20RNA%20sequencing%20in%20the%20prediction%20of%20functional%20and%20structural%20cardiotoxicants%20using%20hiPSC%20card.pdf" "Rosell-Hidalgo_2023_Table1.pdf"
echo "done: $OUT"
