#!/usr/bin/env bash
# download.sh — fetch the GSE247541 GEO submission into source/ (ADR 0030).
#
# method: dataset — NCBI GEO public FTP; no authenticated or browser step.
# source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247541
# The two file names are the ones actually present in source/. This series is a
# provenance-only item: it is outside the Stage 2 CM union scope.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ITEM_DIR="$(dirname -- "$HERE")"
ROOT="$ITEM_DIR"
while [ ! -e "$ROOT/.git" ]; do
  [ "$ROOT" = "/" ] && { echo "download: cannot find the repository root" >&2; exit 2; }
  ROOT="$(dirname -- "$ROOT")"
done
GSE=GSE247541
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/${GSE}"
OUT="$ITEM_DIR/source"
mkdir -p "$OUT"

fetch() { # fetch <url> <name>
  local url="$1" name="$2"
  if [ -s "$OUT/$name" ]; then echo "have  $name"; return 0; fi
  echo "fetch $name"
  curl -fL --retry 3 --retry-delay 2 -o "$OUT/$name.part" "$url"
  mv -- "$OUT/$name.part" "$OUT/$name"
}

# raw counts: GEO supplementary file (.../suppl/)
fetch "${BASE}/suppl/${GSE}_raw_counts.csv.gz" "${GSE}_raw_counts.csv.gz"
# series matrix: GEO matrix file (.../matrix/)
fetch "${BASE}/matrix/${GSE}_series_matrix.txt.gz" "${GSE}_series_matrix.txt.gz"
echo "done: $OUT"
