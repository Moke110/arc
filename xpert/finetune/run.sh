#!/bin/bash
set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate xpert

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

PRETRAINED_MODEL="${1:?Usage: run.sh <pretrained_model_path> [dataset] [split_method] [config] [pretrain_dataset]}"
DATASET="${2:-cdsdb}"
SPLIT_METHOD="${3:-random}"
CONFIG="${4:-finetune/config/cdsdb}"
PRETRAIN_DATASET="${5:-}"

# training hyperparameters
TRAIN_LR="${TRAIN_LR:-0.004}"
WEIGHT_DECAY="${WEIGHT_DECAY:-0.00001}"
BATCH_SIZE="${BATCH_SIZE:-32}"
NUM_EPOCHS="${NUM_EPOCHS:-2}"
INIT_EPOCH="${INIT_EPOCH:-0}"
LOSS_WEIGHT="${LOSS_WEIGHT:-0.2,0.003,0.2,1}"
PATIENCE="${PATIENCE:-50}"
BANDWIDTH="${BANDWIDTH:-0.3}"

# finetune hyperparameters
FINETUNE_LR="${FINETUNE_LR:-0.002}"
FINETUNE_WEIGHT_DECAY="${FINETUNE_WEIGHT_DECAY:-0.00001}"
FINETUNE_PARAMS="${FINETUNE_PARAMS:-all}"
EXCLUDE_KEYS="${EXCLUDE_KEYS:-}"

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
    --pretrain_dataset "${PRETRAIN_DATASET}" \
    --use_gradscaler True \
    --include_cell_idx True \
    --pretrained_model "${PRETRAINED_MODEL}" \
    --train_lr "${TRAIN_LR}" \
    --weight_decay "${WEIGHT_DECAY}" \
    --batch_size "${BATCH_SIZE}" \
    --num_epochs "${NUM_EPOCHS}" \
    --init_epoch "${INIT_EPOCH}" \
    --loss_weight "${LOSS_WEIGHT}" \
    --patience "${PATIENCE}" \
    --bandwidth "${BANDWIDTH}" \
    --finetune_lr "${FINETUNE_LR}" \
    --finetune_weight_decay "${FINETUNE_WEIGHT_DECAY}" \
    --finetune_params "${FINETUNE_PARAMS}" \
    --exclude_keys "${EXCLUDE_KEYS}"
