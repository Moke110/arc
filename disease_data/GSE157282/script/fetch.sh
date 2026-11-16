#!/usr/bin/env bash
set -euo pipefail
ITEM="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ITEM/source"
if [ ! -s "$ITEM/source/GSE157282_RAW.tar" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157282/suppl/GSE157282_RAW.tar" -o "$ITEM/source/GSE157282_RAW.tar.part"
 mv "$ITEM/source/GSE157282_RAW.tar.part" "$ITEM/source/GSE157282_RAW.tar"
fi
if [ ! -s "$ITEM/source/GSE157282_series_matrix.txt.gz" ]; then
 curl -fL --retry 3 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157282/matrix/GSE157282_series_matrix.txt.gz" -o "$ITEM/source/GSE157282_series_matrix.txt.gz.part"
 mv "$ITEM/source/GSE157282_series_matrix.txt.gz.part" "$ITEM/source/GSE157282_series_matrix.txt.gz"
fi
