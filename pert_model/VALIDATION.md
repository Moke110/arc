# Validation log — new model candidates (2026-10-19)

> **Historical log (pre-restructure).** Preprint models (X-Pert, DORA, BaiZe,
> PerturbLDM) were removed from `pert_model/` on 2026-10-19 and are not part of
> Stage 2 setup. The current deployed smoke set is under `stage2/01_setup/`.

This is **setup / smoke validation**, not a benchmark. No model was scored on
CM data; the goal was to prove that the vendored code and released checkpoints
are usable, and to record upstream blockers.

Environment: host Python 3.13 / torch 2.12.0+cu130 for most PyTorch checks;
CMonge uses a Python 3.11 `.venv`, scVIDR a Python 3.8 environment with
scvi-tools 0.13.0, and X-Pert a Python 3.13 venv with two local import shims.
No upstream file was edited inside `pert_model/`.

## Summary

| model | check | result |
|---|---|---|
| PerturbNet | released LINCS + sci-Plex checkpoint load into the in-repo architecture | ✅ all keys matched |
| CMonge | official demo train + evaluate on the bundled sciplex dummy data | ✅ R² 0.886–0.920 |
| Squidiff | 500-gene, drug-structure-conditioned forward pass | ✅ 43.3M params, output `(4, 500)` |
| **scVIDR** | 2,000-gene synthetic train (2 epochs) + continuous-dose prediction | ✅ predictions for doses 1/10/100 |
| **X-Pert** | 2,000-gene drug-mode forward + DoseFiLM gradient | ✅ dose path connected |
| **PerturbLDM** | 2,000-gene VAE encode/decode + conditional diffusion forward | ✅ both shapes `(4, 2000)` |
| **DORA** | 2,000-gene forward with two dose schedules | ✅ dose effect `1.7e-2` |
| **DECANT** | 2,000-gene forward with two dose values | ✅ dose effect `1.6e-2` |
| **BaiZe** | 2,000-gene drug-mode forward + sampling, dose 1 vs 100 | ✅ dose effect `0.202` |
| WAVE | 1-epoch end-to-end train/test on synthetic 978-gene h5ad | ✅ smoke; ❌ no dose gate |
| ExPO | upstream test suite after re-encoding one cp1252 byte | ✅ tests; ❌ 978-only gate |

## PerturbNet

Assets (from the third-party reproduction mirror
`huggingface.co/cyclopeta/PerturbNet_reproduce`, downloaded through
`hf-mirror.com`):

| file | shape / role |
|---|---|
| `assets/example_data_lincs_example.h5ad` | 30,568 cells × 978 genes |
| `assets/example_data_sciplex_example.h5ad` | 31,200 cells × 5,087 genes |
| `assets/models_lincs_..._model_params.pt` | LINCS cINN + ChemicalVAE, 816 tensors |
| `assets/models_sciplex_..._model_params.pt` | sci-Plex cINN + ChemicalVAE |
| `assets/models_sciplex_..._z10_scvi_model_model_params.pt` | sci-Plex cell VAE |
| `assets/models_lincs_..._cellvae_...` | LINCS TF1 cell VAE |
| `assets/pretrained_model_chemicalVAE_model_params_525.pt` | ZINC ChemicalVAE, 4.03M params |

Checkpoint check (PyTorch only; `torchvision` / `scvi` stubbed at import time
because the stable package expects Python 3.7 + TF 1.15):

```
LINCS  cINN: ConditionalFlatCouplingFlow(in=10, cond=196, hidden=1024,
             depth=2, flows=20, bn=True)
             -> 169,843,242 params; `flow.*` all keys matched
             ChemicalVAE(z_dim=196) -> all 53 keys matched
sci-Plex cINN: cond=200 (196 ChemicalVAE + 4 covariates)
             -> `flow.*` all keys matched
             ChemicalVAE(z_dim=196) -> all keys matched
```

Finding: the released weights are structurally intact and match the in-repo
architecture; full inference still needs the Python 3.7/TF 1.15 environment (or
the 0.0.3b1 Python 3.10 branch). The 5,087-gene sci-Plex checkpoint is the
relevant one for Stage 2.

Dose check: the sci-Plex cINN uses `conditioning_dim=200` = 196 ChemicalVAE
dimensions + 3 cell-type one-hots + **one raw numeric dose**; the tutorial builds
`covariates_embedding = concat(cell_type_onehot, adata.obs[['dose']])` and feeds
it to `sample_data`. Because the dose is a raw scalar (not a one-hot), an
arbitrary dose can be supplied at inference.

## CMonge

Created `.venv` with `uv` (Python 3.11, jax 0.4.38, ott-jax 0.5.2). The
official `scripts/demo_train.py` fails immediately because
`configs/demo_config.yml` contains the author's macOS path
`/Users/alicedriessen/Projects/conditional-monge`. After replacing that path
with a local run directory in a **copy** under `/tmp`, the demo completed
end-to-end on CPU:

```
AE train -> CMonge train (100 iters, ~28 s) -> evaluate
abexinostat-10: mean_r2 0.8865
belinostat-10:  mean_r2 0.9199
dacinostat-10:  mean_r2 0.8919
```

This confirms the conditional OT path, RDKit structure embedding, and dose
conditioning all run. The hardcoded path is an upstream packaging bug, not a
model problem.

## Squidiff

```
create_model(gene_size=500, num_layers=2, output_dim=500,
             use_encoder=True, use_drug_structure=True, drug_dimension=1024)
forward(x=(4,500), timesteps=(4,), x_start=(4,500),
        drug_dose=(4,1024), control_feature=(4,500)) -> (4,500)
```

43,342,384 parameters. The drug input is an rFCFP fingerprint multiplied by
`log10(dose+1)`; there is no separate dose covariate. Pretrained checkpoints are
on FigShare (manual).

## WAVE

Synthetic h5ad: 96 train / 32 val / 32 test cells, 978 genes, `layers['unpert_expr']`
and `obs['smiles']`, three valid SMILES.

```
python wave/train.py --outdir . --train_dataset ... --val_dataset ... \
  --test_dataset ... --epochs 1 --batch_size 32 --device cpu
# Test Results | Pearson: 0.0003 | R²: -0.0142 | ΔGene Pearson: 0.7033
# test_pred.h5ad: layers['pred'] (32, 978)
```

Upstream bug: `train_model` saves `best_model.pth` to the **current working
directory**, while `test_model` looks for it under `--outdir`. Running with
`--outdir .` makes the smoke pass; a real run needs a one-line fix. WAVE's model
class hardcodes 978 genes; `GeneVAE(input_dim=3000)` works but `WAVE()` does not
expose the gene width.

## ExPO

`import expo` fails on a clean Python 3 because
`expo/models/expo_model.py` contains a cp1252 en-dash (byte `0x96`) at line 19.
After converting that file to UTF-8 in a scratch copy:

```
python -m pytest tests/test_expo.py -q
# 5 passed, 6 warnings in 43.34s
```

The tests cover imports, config, synthetic data, model instantiation, and dataset
creation at 20/100 genes. No pretrained weights are released.

## scVIDR

Dedicated Python 3.8 environment with scvi-tools 0.13.0, PyTorch 2.4.1 (CPU
wheel was replaced by a GPU-capable build during dependency resolution),
scanpy 1.9.8, geomloss 0.2.5. Synthetic AnnData: 3 cell types × 4 doses
(0/1/10/100) × 20 cells × **2,000 genes**, with the internal categorical
`dose` column required by the CLI.

```
prepare_cont_data(..., 'celltype', 'dose', 'Dose', test_celltype='C',
                  control_dose=0.0)
VIDR(train_adata, linear_decoder=False).train(max_epochs=2, batch_size=64,
                                             early_stopping=False)
model.predict(ctrl_key='0', treat_key='100', cell_type_to_predict='C',
              continuous=True, doses=[0, 1, 10, 100], low_dose=False)
# -> {1: (15, 2000), 10: (15, 2000), 100: (15, 2000)}
```

Training and continuous-dose prediction both work. Caveat: the continuous path
is a **log-linear latent interpolation**
(`delta * log1p(d)/log1p(max(doses))`), not a learned dose-response curve, and
the shipped CLI's `check_dose()` refuses doses absent from `obs`. Treat scVIDR
as a dose-interpolation baseline, not a mechanistic dose model.

## X-Pert (Chen-Li-17, distinct from Stage 1 `xpert/`)

Dedicated Python 3.13 venv with `scperturb==0.1.0`, `einops_exts`, and local
shims for `torchtext.vocab`/`scib` (the repo imports them at module load even
for non-text/non-metrics paths). Drug mode, `dosage_mode_type=4` (DoseFiLM), 2,000
genes:

```
model = TransformerGenerator(ntoken=2001, d_model=32, nhead=4, d_hid=64,
    nlayers=1, nlayers_cls=1, n_cls=1, vocab=..., dropout=0.0,
    pert_mode='drug', output_mode='mlp', gpt_emb_dim=64,
    drug_embed_mode='gpt', dosage_mode_type=4, pert_flag_mode=True,
    use_scgpt_layer=False, feature_dim=2000)
out = model(src, values, flags, mask, pert_emb, pert_mask, doses)
# -> mlm_output (1, 1, 2000)
loss.backward()
# -> DoseFiLM.mlp[-1].weight.grad.norm() = 0.0259
```

The dose path is a real continuous `log10(1+dose)` FiLM modulation. At random
initialisation DoseFiLM is deliberately zero-initialised (γ=1, β=0), so an
untrained model is dose-invariant; the gradient proves the path is connected.
Caveat: the package is alpha-quality — `xpert/__init__.py` exports only the data
loader, `analysis/analyzer.py` references an unexported `PerturbationModel`, and
the chemical workflow lives in notebooks.

## PerturbLDM

No installation needed for the shape check; imported the two model modules
directly. At **2,000 genes**:

```
CellEncoderWithLogvar(latent_dim=32, input_dim=2000, ...)
  encode(x).latents -> (4, 32); decode(...).recon_mu -> (4, 2000)

DenoisingModelDrugDose(latent_dim=32, drug_vocab_size=10, ctrl_dim=2000, ...)
  forward(latents, timesteps, ctrl_expr, drug_id, dose_uM_level)
  -> predict_output (2, 32)
```

Dose enters through `dose_uM_level = log10(dose_uM/0.05)` and an MLP projection;
at random initialisation the dose effect is small (`~1.2e-7`) but the term is
connected. The VAE is the gene-space component; the diffusion model works in
the VAE latent space. No checkpoints or raw data are published.

## DORA

BSD-3, imported `scripts/dora/src/model.py` directly. At **2,000 genes**:

```
DORA(num_genes=2000, num_drugs=5, num_cells=3, genes=(3,2000), hparam=...)
predict(drugs, drugs_feature, cell_feature, indices, cells, drugs_dose)
  with dose schedule [0,1,10] vs [0,10,100]
  -> final gene reconstruction (2, 2000); dose effect L1 = 1.67e-2
```

DORA's dose module is a per-drug `GeneralizedSigmoid` over `log1p(dose)`, applied
as cumulative dose steps. This is the most pharmacological dose parameterisation
in the set. Caveat: no released checkpoint, and the exact gene width depends on
the SciPlex3 preprocessing (the LINCS branch is 978).

## DECANT

Earlier import-only check was too shallow. With the in-repo
`Phase1MechanismTopologyModel` instantiated at **2,000 genes**:

```
model = Phase1MechanismTopologyModel(n_genes=2000, n_programs=50, ...)
forward_view(... dose_log1p=log1p([1,1,1,1]) ...) vs log1p([100,...])
# pred_gene (4, 2000); dose effect L1 = 1.60e-2
```

So the code does accept a continuous dose covariate and a configurable gene
dimension. Caveats: it needs matched treated/control cell bags, dose is treated
as context to suppress rather than a target response axis, and the repository
has no processed data or checkpoints.

## BaiZe

No licence file (public code, all rights reserved — internal evaluation only).
Imported `models/scerso_diffusion.PerturbationDiffusionPredictor` directly at
**2,000 genes**, drug mode, with a 64-d drug feature:

```
model = PerturbationDiffusionPredictor(n_genes=2000, n_perturbations=10,
    drug_dim=64, task_mode='drug', perturb_dim=32, hidden_dims=(64,64),
    timesteps=10, dose_dim=16)
loss, details = model(ctrl, pert, target_rna=target, dose=[1,1],
                      drug_feat=drug, t=[1,1], return_details=True)
pred = model.predict_single(ctrl, pert, dose=[1,1], drug_feat=drug,
                            sample_steps=2)
# pred (2, 2000); dose 1 vs 100 L1 = 0.202
```

The model's dose path is a continuous `Linear(1, dose_dim)` projection and is
clearly active. No main checkpoint or dataset is bundled.

## Reproduce

The exact commands above are intentionally kept outside the vendored trees:
WAVE/ExPO/PerturbNet/PerturbLDM/DORA scratch files are under `/tmp/opencode/`,
CMonge's demo was run in a path-patched copy under `/tmp/opencode/cmonge_run/`,
scVIDR uses `/tmp/opencode/scvidr_env`, and X-Pert uses `/tmp/opencode/xpert_env`
with local `torchtext`/`scib` shims. Re-running the checks does not require
editing any `pert_model/<model>/` source file.
