#!/bin/bash
set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate arc

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

DATASET="${1:-l1000_mdmt_full}"
SPLIT_METHOD="${2:-pretrain_split}"
# Pretrain config + base weight live in the consuming workflow (ADR 0008);
# the L1000 pretrain h5ads live in pert_data/l1000-xpert/ (dataset home).
CONFIG="${3:-../../workflow/xpert-fda-miri-v1/03_xpert_pretrain/config/l1000}"

# training hyperparameters
TRAIN_LR="${TRAIN_LR:-0.004}"
WEIGHT_DECAY="${WEIGHT_DECAY:-0.00001}"
BATCH_SIZE="${BATCH_SIZE:-64}"
NUM_EPOCHS="${NUM_EPOCHS:-2500}"
INIT_EPOCH="${INIT_EPOCH:-70}"
LOSS_WEIGHT="${LOSS_WEIGHT:-0.2,0.003,0.2,1}"
PATIENCE="${PATIENCE:-50}"
BANDWIDTH="${BANDWIDTH:-0.3}"

# drug embeddings to use (comma-separated: unimol,KPGT,morgan,smi,HG)
DRUG_EMBEDS="${DRUG_EMBEDS:-unimol,HG}"
# primary molecular feature
DRUG_FEAT="${DRUG_FEAT:-unimol}"

python -m src.train \
    --mode train \
    --model XPert \
    --config "${CONFIG}" \
    --drug_feat "${DRUG_FEAT}" \
    --drug_embeds "${DRUG_EMBEDS}" \
    --split_method "${SPLIT_METHOD}" \
    --dataset "${DATASET}" \
    --use_gradscaler True \
    --include_cell_idx True \
    --train_lr "${TRAIN_LR}" \
    --weight_decay "${WEIGHT_DECAY}" \
    --batch_size "${BATCH_SIZE}" \
    --num_epochs "${NUM_EPOCHS}" \
    --init_epoch "${INIT_EPOCH}" \
    --loss_weight "${LOSS_WEIGHT}" \
    --patience "${PATIENCE}" \
    --bandwidth "${BANDWIDTH}"
