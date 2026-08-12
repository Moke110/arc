#!/usr/bin/env bash
set -euo pipefail

./.venv/bin/python -m src.train_transigen \
  --training_data_root training/data \
  --dataset_name GSE92742 \
  --pretrained_vae1_path training/pretrained/vae1.pt \
  --pretrained_vae2_path training/pretrained/vae2.pt \
  --use_pretrained_vae true \
  --output_root training/models \
  --dev cuda:2 \
  --n_epochs 500 \
  --n_latent 100 \
  --molecule_feature_embed_dim 400 \
  --batch_size 64 \
  --learning_rate 1e-3 \
  --beta 0.1 \
  --dropout 0.1 \
  --weight_decay 1e-5 \
  --patient_cons 10 \
  --patience_turn 30 \
  --train_flag true \
  --eval_metric false \
  --predict_profile false
