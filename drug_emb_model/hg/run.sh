#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python -m src.pretrain_hg \
    --include_edges PPI DTI DDS \
    --hg_path "${SCRIPT_DIR}/data/" \
    --num_layers 3 \
    --hidden_dim 256 \
    --out_dim 256 \
    --num_neighbors 10 5 3 \
    --batch_size 1 \
    --epochs 300 \
    --patience 30 \
    --conseq_patience 3 \
    --turn_patience 10 \
    --num_neg_samples 3 \
    --temperature 0.5 \
    --device cuda:0 \
    "$@"
