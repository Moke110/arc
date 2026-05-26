#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 在这里固定要转换的 CSV 路径（相对仓库根目录）。
SMILES_CSV="feature_extract/INO/smiles.csv"

python "$SCRIPT_DIR/src/smiles_to_kpgt_pickle.py" \--smiles_csv "$SMILES_CSV" "$@"
