#!/usr/bin/env bash
# fetch.sh — download the GSE141910 GEO supplementary files into source/ (ADR 0030).
#
# method: dataset — NCBI GEO public FTP (HTTPS mirror); no authenticated or
# browser step. Exactly two upstream files are acquired, byte-unmodified:
#   GSE141910_RAW.tar   60,026,880 B  366 per-sample ENSG count CSVs
#   filelist.txt            22,419 B  upstream listing of those 366 CSVs
# Both live under .../suppl/. The series matrix is NOT downloaded: the
# case/control labels are the two group tokens in the per-sample file names
# (the same titles the series matrix carries) and they are pinned in
# script/build_deg.R; see author_log.json "notes".
#
# Usage:  cd disease_data/GSE141910 && bash script/fetch.sh
# Then:   bash script/verify.sh   (offline, sha256sum -c ../manifest.sha256)
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ITEM_DIR="$(dirname -- "$HERE")"
ROOT="$ITEM_DIR"
while [ ! -e "$ROOT/.git" ]; do
  [ "$ROOT" = "/" ] && { echo "fetch: cannot find the repository root" >&2; exit 2; }
  ROOT="$(dirname -- "$ROOT")"
done
GSE=GSE141910
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/${GSE}"
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

fetch "${BASE}/suppl/${GSE}_RAW.tar" "${GSE}_RAW.tar" 60026880
fetch "${BASE}/suppl/filelist.txt" "filelist.txt" 22419
echo "done: $OUT"
