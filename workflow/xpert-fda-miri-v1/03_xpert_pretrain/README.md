# Step 03 — XPert pretrain (reference)

This step owns the **author-given XPert pretrain unit** used as the base weight
for the DToxS finetune (step 05) and the DTOXS-split benchmark (step 04):

- `pretrain_mdmt_full_50_epoch.pth` — the used pretrained XPert checkpoint
  (L1000, 50 epochs). This is the exact weight the workflow's finetune loads by
  default.
- `config/l1000.yaml` — the pretrain recipe that produced the weight. Data
  roots are repointed after the workspace refactor (ADR 0008): the three L1000
  pretrain h5ads now live in the dataset home `pert_data/l1000/data/`, and the
  drug embeddings (UniMol array, SMILES mapping, author HG reference) are
  workflow-owned in `input/`.

The two unused pretrain checkpoints (`pretrain_mdmt_full_200_epoch.pth`,
`pretrain_mdmt_new.pth`) were deleted during the refactor; only the used weight
is preserved.

`pert_model/xpert/pretrain/run.sh` remains in the model library as the
pretrain recipe; its default config now resolves here
(`../workflow/xpert-fda-miri-v1/03_xpert_pretrain/config/l1000`).

Note: pretraining is not re-run during screening — steps 04/05 consume the
weight from this step directory.
