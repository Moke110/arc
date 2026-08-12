#!/usr/bin/env bash
#  --seed 364039 \

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

"${ROOT_DIR}/.venv/bin/python" -m src.prediction \
  --model_path training/models/LINCS2020_260422-1953/best_model.pt \
  --predict_data_root prediction/INO/ \
  --compound_info screening_compound.csv \
  --compound_emb KPGT_emb.pickle \
  --cell HUVEC \
  --seed 364039 \
  --dev cuda:0 \
  --batch_size 64 \
  --num_workers 2
