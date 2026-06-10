#!/bin/bash
# XPert drug perturbation prediction
#
# Predicts transcriptional perturbation (deg = X_pred - X_ctl) for a list of
# drugs given a baseline cell transcriptome and model weights.
#
# Directory layout under prediction_root:
#   {prediction_root}/
#   ├── process/           ← intermediate files (UniMol cache, skipped drugs)
#   ├── {ctl_csv}          ← input: baseline transcriptome (978 landmark genes)
#   ├── {drug_csv}         ← input: drug list (must have 'canonical_smiles' column)
#   └── {output_csv}       ← output: predicted deg matrix
#
# Usage:
#   bash predict/predict.sh \
#       --prediction_root predict/fda-cm \
#       --ctl_csv cm_ctl_landmark.csv \
#       --drug_csv fda_drugs.csv \
#       --output_csv deg-perturb.csv \
#       --model_path finetune/output/dtoxs_warm/ \
#       --dose_idx 0 --time_idx 2

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
PREDICTION_ROOT=""
CTL_CSV=""
DRUG_CSV=""
OUTPUT_CSV="deg-perturb.csv"
DOSE_IDX=0
TIME_IDX=2
MODEL_PATH=""
DEVICE="cuda:0"
BATCH_SIZE=64

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prediction_root) PREDICTION_ROOT="$2"; shift 2 ;;
        --ctl_csv)         CTL_CSV="$2";         shift 2 ;;
        --drug_csv)        DRUG_CSV="$2";        shift 2 ;;
        --output_csv)      OUTPUT_CSV="$2";      shift 2 ;;
        --dose_idx)        DOSE_IDX="$2";        shift 2 ;;
        --time_idx)        TIME_IDX="$2";        shift 2 ;;
        --model_path)      MODEL_PATH="$2";      shift 2 ;;
        --device)          DEVICE="$2";          shift 2 ;;
        --batch_size)      BATCH_SIZE="$2";      shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate required args
: ${PREDICTION_ROOT:?Missing --prediction_root}
: ${CTL_CSV:?Missing --ctl_csv}
: ${DRUG_CSV:?Missing --drug_csv}
: ${MODEL_PATH:?Missing --model_path}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
cd "$PROJECT_ROOT"

python -m src.predict \
    --prediction_root "$PREDICTION_ROOT" \
    --ctl_csv "$CTL_CSV" \
    --drug_csv "$DRUG_CSV" \
    --output_csv "$OUTPUT_CSV" \
    --dose_idx "$DOSE_IDX" \
    --time_idx "$TIME_IDX" \
    --model_path "$MODEL_PATH" \
    --device "$DEVICE" \
    --batch_size "$BATCH_SIZE"
