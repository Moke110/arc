# Download — PerturbNet

**method: repo + HF-mirror asset** — upstream clone (GPL-3.0) plus selected files
from a third-party reproduction mirror.

## Repository

```bash
git clone https://github.com/welch-lab/PerturbNet pert_model/PerturbNet/source
git -C pert_model/PerturbNet/source checkout e1375b2ca71c87d5b7277c45648ed06fe335a38d
# vendored read-only: the nested .git directory was removed on 2026-10-19
```

## Assets (not author-released)

The released checkpoints do not exist in the upstream repository; the Stage 2
assets are selected files from the reproduction mirror
<https://huggingface.co/cyclopeta/PerturbNet_reproduce>, fetched through the
`hf-mirror.com` endpoint and kept at `pert_model/PerturbNet/source/assets/`:

| file | role |
|---|---|
| `example_data_lincs_example.h5ad` | 30,568 cells x 978 genes |
| `example_data_sciplex_example.h5ad` | 31,200 cells x 5,087 genes |
| `models_lincs_*_cinn_*_model_params.pt` | LINCS cINN (664 MB) |
| `models_sciplex_*_cinn_*_model_params.pt` | sci-Plex cINN (664 MB) |
| `models_sciplex_*_z10_scvi_model_model_params.pt` | sci-Plex cell-VAE (13 MB) |
| `models_lincs_*_cellvae_*` | LINCS TF1 cell-VAE (17 MB) |
| `pretrained_model_chemicalVAE_model_params_525.pt` | ZINC ChemicalVAE (16 MB) |

`source/ASSETS.md` (shipped by the mirror) holds the per-file sizes and SHA256s.

## Notes

Deployed model; requires Python 3.7 + TensorFlow 1.15 for the stable branch (or
the 0.0.3b1 py3.10 branch without TF). Run card `stage2/01_setup/perturbnet/`.
