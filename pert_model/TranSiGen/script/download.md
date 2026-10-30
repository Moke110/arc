# Download — TranSiGen

**method: repo + results asset** — a nested git clone (MIT) plus a separately
obtained pretrain-checkpoint/result set.

## Repository

The upstream clone was kept as a nested git repository and that `.git` directory
was removed on 2026-10-30 after recording its commit:

```
repo:   https://github.com/myzhengSIMM/TranSiGen
branch: main
commit: e7582e7794f13d4dfddfb15e95e3ba7cd24bc07a
date:   2025-01-21  (last commit "Add files via upload")
```

To reproduce:

```bash
git clone --branch main https://github.com/myzhengSIMM/TranSiGen \
  pert_model/TranSiGen/source
git -C pert_model/TranSiGen/source checkout e7582e7794f13d4dfddfb15e95e3ba7cd24bc07a
cd pert_model/TranSiGen/source && rm -rf .git
```

`source/src/`, `source/data/` (LINCS2020 splits/embeddings, PRISM,
`TranSiGen.jpg`), `LICENSE` and `README.md` are clone content.

## Results asset

`source/results/` holds the pretrain checkpoints and the published prediction
outputs, materialized from the external weight archive on 2026-10-30 (ticket 02):

- `trained_model_shRNA_vae_x1/best_model.pt` (9,191,951 B)
- `trained_model_shRNA_vae_x2/best_model.pt` (9,191,741 B)
- `trained_models_164_cell_smiles_split/364039/feature_KPGT_init_pretrain_shRNA/best_model.pt` (22,441,478 B)
- `5.Drug_response_prediction/`, `6.Phenotype_based_drug_repurposing/` result tables

## Notes

Reference-only, kept for Stage 1 provenance; not a Stage 2 entrant.
