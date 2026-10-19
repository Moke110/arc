# Stage 2 model assets

The former repository-root `download/` directory was integrated here on
2026-10-19. The h5ad files are the CPA Drive datasets; the zips are the original
weight archives. All are gitignored; this manifest is the tracked provenance
record.

## Datasets

| file | size | sha256 | consumer |
|---|---:|---|---|
| `cpa/data/kang_normalized_hvg.h5ad` | 544,912,376 | `f04d4d31cc5ff56588e34cb305cdf52b3eddaad273192857a95dd7f5bea1c437` | CPA smoke |
| `cpa/data/Norman2019_normalized_hvg.h5ad` | 837,969,230 | `70f9f81cff6621b06f98fe150fd8e7e639fbd252b8d4fe1ad98450d54a6bb917` | CPA (future; Norman checkpoint) |
| `biolord/data/GSM_new.h5ad` | 178,640,710 | `11b4363d46c00bcb9109c232b4f18411d01889d7a218f459a7d3feaa9b8f23da` | biolord smoke |
| `chemcpa/data/combo_sciplex_prep_hvg_filtered.h5ad` | 490,412,774 | `347f3a0dbfe73fa143adcf0eee54dc21a3a97129e0a77d946800cb8981b32141` | chemCPA (future retraining; not used by the load-only smoke) |

## Setup-stage assets (2026-10-19)

- `perturbnet/assets/` — selected LINCS/sci-Plex checkpoints and examples; see
  `perturbnet/ASSETS.md` for sizes/SHA256.
- `scvidr/source/` — optional `data.zip` (manual; see `../download/manual.md`).
- `squidiff/weights/` — optional FigShare checkpoints (manual).
- `cmonge/data/` — bundled demo data (no download).

## Source archives

| file | size | sha256 |
|---|---:|---|
| `cpa/source/Pretrained_Models-20261018T124453Z-1-001.zip` | 59,157,883 | `ab0d6ef9cbd000aa1a150c0475671cc2eb1f1c3ab8b67aa1d9fb2145f896967d` |
| `chemcpa/source/chemCPA_models.zip` | 254,487,527 | `322606b8a209de44228cf51d9930e4173b21a79dfc780fab3c1a2f1e916bec8f` |

## Extracted weights

- `cpa/weights/Pretrained_Models/{Kang,Norman2019}/`
- `chemcpa/weights/chemCPA_models/*.pt` (12 checkpoints; L1000 977/2,000 genes)
- `prnet/checkpoint/{lincs,sciplex}_best_epoch_all.pt` (shipped with the repo)

## Provenance

- CPA/pretrained data: <https://drive.google.com/drive/folders/1yFB0gBr72_KLLp1asojxTgTqgz6cwpju>
- chemCPA weights: <https://f003.backblazeb2.com/file/chemCPA-models/chemCPA_models.zip>
- PRnet: <https://github.com/Perturbation-Response-Prediction/PRnet>
- biolord: <https://github.com/nitzanlab/biolord>

Per-model `VENDORED.txt` files record repository commits and licenses.
