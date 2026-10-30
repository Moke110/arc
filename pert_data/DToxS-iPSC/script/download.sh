#!/usr/bin/env bash
# download.sh — fetch the DToxS iPSC-CM GEO deposit into source/ (ADR 0030).
#
# method: dataset — NCBI GEO public FTP; no authenticated or browser step.
# source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217424
# GSE217424 is the SuperSeries of two sub-series:
#   GSE217421 — Conv arm (conventional/hiPSC-CM RNA-seq)
#   GSE217423 — DGE arm (3'-DGE RNA-seq)
# 17 files are actually present in source/ (9 + 8); :Zone.Identifier sidecars are
# not downloaded.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ITEM_DIR="$(dirname -- "$HERE")"
ROOT="$ITEM_DIR"
while [ ! -e "$ROOT/.git" ]; do
  [ "$ROOT" = "/" ] && { echo "download: cannot find the repository root" >&2; exit 2; }
  ROOT="$(dirname -- "$ROOT")"
done
OUT="$ITEM_DIR/source"
mkdir -p "$OUT"

fetch() { # fetch <url> <name>
  local url="$1" name="$2"
  if [ -s "$OUT/$name" ]; then echo "have  $name"; return 0; fi
  echo "fetch $name"
  curl -fL --retry 3 --retry-delay 2 -o "$OUT/$name.part" "$url"
  mv -- "$OUT/$name.part" "$OUT/$name"
}

# GSE217421 (Conv arm) — 9 files
GSE=GSE217421
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/${GSE}"
for f in Conv-DEG-Calc-FDR-0.1.tar.gz Conv-RNAseq-Configs.tsv.gz \
         Conv-RNAseq-Experiment-Design.tsv.gz Conv-RNAseq-Params.tar.gz \
         Conv-RNAseq-Read-Counts.tsv.gz Program-Source-Codes.tar.gz \
         README.pdf.gz Reference-Library.tar.gz; do
  fetch "${BASE}/suppl/${GSE}_${f}" "${GSE}_${f}"
done
fetch "${BASE}/matrix/${GSE}_series_matrix.txt.gz" "${GSE}_series_matrix.txt.gz"

# GSE217423 (DGE arm) — 8 files
GSE=GSE217423
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/${GSE}"
for f in DGE-DEG-Calc-FDR-0.1.tar.gz DGE-RNAseq-Configs.tsv.gz \
         DGE-RNAseq-Experiment-Design.tsv.gz DGE-RNAseq-Params.tar.gz \
         DGE-RNAseq-Read-Counts.tsv.gz Program-Source-Codes.tar.gz \
         Reference-Library.tar.gz; do
  fetch "${BASE}/suppl/${GSE}_${f}" "${GSE}_${f}"
done
fetch "${BASE}/matrix/${GSE}_series_matrix.txt.gz" "${GSE}_series_matrix.txt.gz"
echo "done: $OUT"
