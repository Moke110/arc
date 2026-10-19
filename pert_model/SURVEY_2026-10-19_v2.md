# Virtual-cell drug-perturbation survey — revised gates (2026-10-19)

This supersedes `SURVEY_2026-10-19.md` for admission decisions. The original
four-gate survey is kept for provenance.

> **Deployment decision (2026-10-19):** preprint models are excluded from Stage 2
> setup. X-Pert, DORA, BaiZe, and PerturbLDM remain in this document as a
> literature record only; their code was removed from `pert_model/` and they have
> no run card under `stage2/01_setup/`.

## Revised gates

| gate | requirement |
|---|---|
| **G1** | gene coverage is **not fixed at 978**; roughly **2,000–20,000** genes (dataset HVG space or configurable width is fine) |
| **G2** | **dose handling / generalisation** — dose is an input or covariate and can be varied at inference |
| **G3** | **public code** — a working public repository; a LICENSE file is recorded but not a hard gate for exploration |
| — | pretrained weights are optional |

## Re-audit of the already-downloaded models

| model | G1 genes | G2 dose | G3 code | verdict now | action |
|---|---|---|---|---|---|
| **PerturbNet** | 978 / **5,087** | cINN covariate, continuous | GPL-3.0 | **PASS** | keep as primary |
| **CMonge** | per-dataset HVG (not 978-fixed) | **categorical `drug-dose` condition** (e.g. `givinostat-10`, `-100`); no continuous dose input | MIT | **BORDERLINE G2** | keep as condition-transfer baseline, not a dose-response model |
| **Squidiff** | **configurable `gene_size`** | dose enters as `rFCFP · log10(dose+1)` | BSD-3 | **PASS (weak dose)** | keep; dose is scaling, not a learned curve |
| **DECANT** | **`n_genes` parameter** | `dose_log1p` in context encoder; dose changes output | MIT | **PASS (contextual dose)** | keep as representation model; no data/weights |
| **scVIDR** | **5,000 HVG** in paper; 2,000 verified | continuous latent interpolation `δ·log1p(d)/log1p(max)` | MIT | **PASS/BORDERLINE** | add; dose is interpolated, CLI blocks unseen doses |
| **X-Pert** | 2,000 verified; `var_num=5000` | continuous `log10(1+dose)` → DoseFiLM | MIT | **PASS/BORDERLINE — preprint, excluded** | not deployed; code removed |
| **WAVE** | 978 default (3,000 possible) | **no dose input** | MIT | **FAIL G2** | demote to reference |
| **ExPO** | **978** (constant default; no 2k+ run) | continuous dose × time | MIT | **FAIL G1** | demote to reference / protocol only |
| **CPA** (Stage 2) | 20,570 VAE dim | dose covariate | BSD-3 | **PASS** | already admitted |
| **chemCPA** (Stage 2) | 20,570 after retraining | dose covariate | MIT | **PASS** | already admitted |
| **biolord** (Stage 2) | 20,570 VAE dim | ordered dose attribute | BSD-3 | **PASS** | already admitted |
| **PRnet** (Stage 2) | configurable x_dimension | native dose encoder | Apache-2.0 | **PASS** | already admitted |

Key correction: **DECANT is not a skeleton-only import failure**. The in-repo
`Phase1MechanismTopologyModel` has a configurable `n_genes`, a `ContextEncoder`
that consumes `dose_log1p`, and a gene decoder. It runs at 2,000 genes and the
output changes with dose. The repository still ships no processed data or
checkpoints, and its design treats dose as context to suppress rather than a
response axis to predict.

A second correction: **CMonge's dose is a categorical condition label**
(`givinostat-10`, `givinostat-100`, …), not a continuous covariate. It handles
multiple observed doses but cannot query an unseen dose without a new condition;
treat it as a condition-transfer baseline, not a dose-response model.

WAVE and ExPO remain in `pert_model/` for comparison, but neither is an entrant
under the revised gates: WAVE has no dose input, ExPO's published scope is 978.

## Newly added candidates

| model | paper | repo | licence | G1 genes | G2 dose | weights/data | verdict |
|---|---|---|---|---|---|---|---|
| **PerturbLDM** | bioRxiv 2026, `10.64898/2026.08.07.743610` | `davidroad/PerturbLDM` | **MIT** | Tahoe-100M panel (tens of thousands); VAE `input_dim` from data | **continuous** `log10(dose_uM/0.05)` → MLP | no checkpoints/raw data | **PASS — preprint, excluded** |
| **DORA** | “Cell Fate Prediction Using Dose-Response Autoencoder Model” (2025) | `LBiophyEvo/DORA` | **BSD-3** | SciPlex3 gene space (LINCS branch 978) | **per-drug generalized sigmoid** over `log1p(dose)`, cumulative steps | no DORA checkpoint | **PASS — preprint, excluded** |
| **scVIDR** | Patterns 2023, `10.1016/j.patter.2023.100817` | `BhattacharyaLab/scVIDR` | **MIT** | 5,000 HVG (2,000 verified) | continuous interpolation; single/multi-dose modes | `data.zip` on Google Drive includes pretrained models | **PASS/BORDERLINE** |
| **X-Pert** | bioRxiv `10.1101/2025.11.13.688367` | `Chen-Li-17/X-Pert` | **MIT** | 2,000 verified; `var_num=5000` | **continuous** `log10(1+dose)` → DoseFiLM | no weights; FigShare data | **PASS/BORDERLINE — preprint, excluded** |
| **DECANT** | Bioinformatics 2026, `10.1093/bioinformatics/btag662` | `bliulab/DECANT` | **MIT** | configurable `n_genes` (2,000 verified) | `dose_log1p` context input | none | **PASS/BORDERLINE** |
| **BaiZe** | bioRxiv 2026, `10.64898/2026.07.15.738608` | `simonzqw/BaiZe` | **no LICENSE** (public code, all rights reserved) | full AnnData transcriptome (2,000 verified) | **continuous** `log1p(dose)` → `Linear(1,32)` | no main checkpoint; no data | **PASS\* — preprint, excluded** |

Why these are interesting:

- **PerturbLDM** is the cleanest new fit: MIT, latent diffusion, continuous dose
  as an independent scalar, and a VAE whose gene width follows the training
  matrix. It is the only new model that combines a wide gene space with a
  genuinely continuous dose covariate and a permissive licence.
- **DORA** has the most pharmacological dose parameterisation: a learnable
  generalized sigmoid per drug, applied over cumulative dose steps. It can be
  queried at arbitrary doses, and the code is BSD-3.
- **scVIDR** is a practical dose-interpolation baseline. Its `predict(
  continuous=True, doses=[...])` returns a prediction per dose, but the curve is
  a log-linear scaling of one latent delta. The CLI refuses doses absent from
  `obs`.
- **X-Pert** has a clean continuous DoseFiLM and a 2,000-gene drug-mode forward.
  The package is early-stage and the chemical tutorials are L1000-derived, so
  measure the drug-branch width yourself before adopting it.
- **DECANT** is a representation model, not a generative dose-response model.
  Dose is an input, but the objective suppresses context; use it for
  mechanism-aligned embeddings, not for dose-response curves.
- **BaiZe** is architecturally the closest to the full requirement: Δx
  prediction, continuous dose, SMILES structure, optional ATAC/time/species, and
  a configurable gene space. The blocker is legal, not technical — there is no
  LICENSE file, so it is internal-evaluation only until the authors add one.

## Checked and excluded

| model | reason |
|---|---|
| **PertFlow** (`wangmengbo/PertFlow`) | Public code but **no LICENSE**, no README; dose/concentration is present in `encoders.py`/`dataloader.py` but the core training path is undocumented. Not validated. |
| **AetherCell** (`Wenyuan-AI4science/AetherCell`) | 978 or 10,085 genes, but **no dose** in the code/API; custom non-OSI licence. Fail G2. |
| **PertDiT** (`wangkekekeke/PertDiT`) | **978** (`d_pre: 978` in every config). Fail G1. |
| **PertDiff** (`Panda-myj/PertDiff`) | **no dose** (0 hits). Fail G2. |
| **scNODE** (`rsinghlab/scNODE`) | 2,000 HVG, MIT, but the continuous axis is **time**, not dose. Fail G2. |
| **trVAE** (`theislab/trVAE`) | dose would be one-hot categorical (`to_categorical`); no continuous covariate. Fail G2. |
| **CellPert** (`QSong-github/CellPert`) | No dose anywhere in `src/`. Fail G2. |
| **scCycleMol** (`Dean0131/scCycleMol`) | Repo empty. Fail G3. |
| **WAVE / ExPO** | see re-audit above. |
| **CellCap / W1OT / scPerb / scPILOT / MultiPert / GPO-VAE / Stable-Shift / AROMA** | no dose or genetic-only; fail G2. |

## Validation performed

All newly admitted code was downloaded into `pert_model/` and exercised at
2,000 genes; details and commands are in `VALIDATION.md`. The preprint rows are
now historical: their code was removed on 2026-10-19 and they are not part of
Stage 2 setup. Summary:

| model | validation |
|---|---|
| PerturbLDM | 2,000-gene VAE encode/decode; conditional diffusion forward; dose path connected |
| DORA | 2,000-gene forward; two dose schedules; dose effect `1.67e-2` |
| scVIDR | 2,000-gene training (2 epochs); continuous predictions for doses 1/10/100 |
| X-Pert | 2,000-gene drug-mode forward; DoseFiLM gradient `2.59e-2` |
| DECANT | 2,000-gene forward; dose effect `1.60e-2` |
| BaiZe | 2,000-gene drug-mode forward; dose 1 vs 100 effect `0.202` |

## Recommendation

Under the revised gates, the working set is:

```
PerturbNet     5,087 genes   continuous dose     GPL-3.0   deployed
scVIDR         5,000 HVG     interpolated dose   MIT       deployed
DECANT         configurable  contextual dose     MIT       deployed
Squidiff       configurable  scaled dose         BSD-3     deployed
CMonge         dataset HVG   categorical dose     MIT       deployed
CPA/chemCPA/biolord/PRnet   already admitted Stage 2
PerturbLDM / DORA / X-Pert / BaiZe   preprints — excluded, code removed
```

Demote WAVE and ExPO to reference/baseline status. Preprint models are excluded
from deployment; PertFlow stays on the watch list (no LICENSE, no README).

The revised gates do **not** remove the central engineering constraint: no model
ships a ready 20,570-gene CM checkpoint. The realistic paths remain (a) retrain a
configurable model on the CM gene space and (b) use a 978/5k model with an
expansion/imputation step. The deployed set adds dose parameterisations
(scVIDR interpolation, PerturbNet continuous covariates) rather than a ready
whole-transcriptome checkpoint.
