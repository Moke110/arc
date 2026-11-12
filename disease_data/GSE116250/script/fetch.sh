#!/usr/bin/env bash
# fetch.sh — download the GSE116250 GEO supplementary file into source/ (ADR 0030).
#
# method: dataset — NCBI GEO public FTP (HTTPS mirror); no authenticated or
# browser step. Exactly one upstream file is acquired, byte-unmodified:
#   GSE116250_rpkm.txt.gz   7,487,245 B
#     gene-level RPKM matrix: `Gene` (ENSG) + `Common_name` (symbol) + 64
#     sample columns. The series matrix is NOT downloaded: the case/control
#     labels are the GSM column identifiers, and the GSM -> group map extracted
#     from the GEO series matrix is pinned in script/build_deg.R; see
#     author_log.json "notes".
#
# Usage:  cd disease_data/GSE116250 && bash script/fetch.sh
# Then:   bash script/verify.sh   (offline, sha256sum -c ../manifest.sha256)
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ITEM_DIR="$(dirname -- "$HERE")"
ROOT="$ITEM_DIR"
while [ ! -e "$ROOT/.git" ]; do
  [ "$ROOT" = "/" ] && { echo "fetch: cannot find the repository root" >&2; exit 2; }
  ROOT="$(dirname -- "$ROOT")"
done
GSE=GSE116250
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/${GSE}"
OUT="$ITEM_DIR/source"
mkdir -p "$OUT"

fetch() { # fetch <url> <name> <expected_bytes>
  local url="$1" name="$2" want="$3"
  if [ -s "$OUT/$name" ]; then echo "have  $name"; return 0; fi
  echo "fetch $name"
  curl -fL --retry 3 --retry-delay 2 -o "$OUT/$name.part" "$url"
  local got
  got="$(stat -c %s "$OUT/$name.part")"
  if [ "$got" != "$want" ]; then
    echo "fetch: $name is $got B, expected $want B" >&2
    rm -f -- "$OUT/$name.part"
    exit 1
  fi
  mv -- "$OUT/$name.part" "$OUT/$name"
}

fetch "${BASE}/suppl/${GSE}_rpkm.txt.gz" "${GSE}_rpkm.txt.gz" 7487245
echo "done: $OUT"
