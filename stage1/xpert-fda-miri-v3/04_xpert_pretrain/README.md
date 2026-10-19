# Step 04 — XPert pretrain (reference)

This step owns the **author-given XPert pretrain unit** used as the base weight
for the union-CM finetune (step 06) and the warm-start benchmark (step 05):

- `pretrain_mdmt_full_50_epoch.pth` — the used pretrained XPert checkpoint
  (L1000, 50 epochs). This is the exact weight the workflow's finetune loads by
  default.
- `config/l1000.yaml` — the pretrain recipe that produced the weight. Data
  roots are repointed after the workspace refactor (ADR 0008): the three L1000
  pretrain h5ads now live in the dataset home `pert_data/l1000-xpert/data/`, and the
  drug embeddings (UniMol array, SMILES mapping, author HG reference) are
  workflow-owned in `input/`.

The two unused pretrain checkpoints (`pretrain_mdmt_full_200_epoch.pth`,
`pretrain_mdmt_new.pth`) were deleted during the refactor; only the used weight
is preserved.

`pert_model/xpert/pretrain/run.sh` remains in the model library as the
pretrain recipe; its default config now resolves here
(`stage1/xpert-fda-miri-v3/04_xpert_pretrain/config/l1000`).

Note: pretraining is not re-run during screening — steps 05/06 consume the
weight from this step directory.
