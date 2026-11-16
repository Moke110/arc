#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p source
if bash script/verify.sh >/dev/null 2>&1; then exit 0; fi
curl --fail --location --continue-at - --retry 4 --connect-timeout 30 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298023/suppl/GSE298023_HH_dataset_seurat_object_2026.rds' --output source/GSE298023_HH_dataset_seurat_object_2026.rds
bash script/verify.sh
