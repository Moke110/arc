# pert_model — vendored drug-perturbation models

The perturbation-prediction model units the pipeline compares, split into one
directory per **upstream source** (项) under the external-data provenance
contract ([ADR 0030](../docs/adr/0030-external-data-provenance-and-step-layout.md)).
Each item is:

```
<Item>/
├── script/          downloader (download.md) + offline verifier (verify.sh)   [tracked]
├── source/          upstream clone + downloaded assets, bytes unchanged       [gitignored]
├── local/           declared local additions (only where they exist)          [gitignored]
└── author_log.json  source URL/commit/licence, asset origins, verification    [tracked]
```

`manifest.sha256` (per-file hashes over `source/` only) exists for every item — 453.0 GB hashed
across the four categories, 27/27 verified. `script/verify.sh` checks it offline, and a
missing manifest is a failure. Run `script/verify_all.sh` from the repository root to check
every item.

## The fourteen historical 项

| 项 | upstream | licence | stage | what it is |
|---|---|---|---|---|
| `PRnet/` | [PRnet](https://github.com/Perturbation-Response-Prediction/PRnet) `f19174bd` | Apache-2.0 | historical setup | PRnet reaction/response predictor; ships its own 978 LINCS + 5,000 Sci-Plex checkpoints. |
| `CPA/` | [CPA](https://github.com/theislab/CPA) `fbd7c025` | BSD-3-Clause | historical setup | Compositional Perturbation Autoencoder; Kang + Norman2019 pretrained models on the authors' Drive. |
| `chemCPA/` | [chemCPA](https://github.com/theislab/chemCPA) `43e830eb` | MIT | historical setup | Chemistry-aware CPA; public 977/2,000-gene L1000 ComPert weights, not retrained on CM data. |
| `biolord/` | [biolord](https://github.com/nitzanlab/biolord) `b7688790` | BSD-3-Clause | historical setup | Disentangled representation model; no public weights (per-dataset training). |
| `PerturbNet/` | [PerturbNet](https://github.com/welch-lab/PerturbNet) `e1375b2c` | GPL-3.0 | historical setup | Conditional INN + ChemicalVAE; 5,087-gene sci-Plex checkpoint with a raw continuous dose covariate. |
| `scVIDR/` | [scVIDR](https://github.com/BhattacharyaLab/scVIDR) `bcb6a96c` | MIT | historical setup | VAE with single-/multi-dose modes; 5,000 HVG, log-linear latent dose interpolation. |
| `Squidiff/` | [Squidiff](https://github.com/siyuh/Squidiff) `abdfc27d` | BSD-3-Clause | historical setup | Diffusion model; configurable gene space, drug input is rFCFP × log10(dose+1). |
| `CMonge/` | [Conditional-Monge](https://github.com/AI4SCR/Conditional-Monge) `195f5cd1` | MIT | historical setup | Conditional Monge-gap transport; RDKit structure + categorical `drug-dose`, bundled demo data. |
| `DECANT/` | [DECANT](https://github.com/bliulab/DECANT) `e50642a9` | MIT | historical setup | Mechanism/context-decoupled drug representation; dose is context to suppress, no shipped data/weights. |
| `WAVE/` | [WAVE](https://github.com/ZJUFanLab/WAVE) `e27e453b` | MIT | reference | GeneDrugFusion; 978 default, **no dose input** → fails revised G2. |
| `ExPO/` | [ExPO](https://github.com/MLBC-lab/ExPO) `7e88cbae` | MIT | reference | Exposure-aware model; **978 only** in the published scope → fails revised G1. |
| `State/` | [state](https://github.com/ArcInstitute/state) `9bbfe78a` | CC BY-NC-SA 4.0 / State Model Non-Commercial | reference | Arc State transition model (SE-600M + ST-HVG-Tahoe); Tahoe one-hot limits cold-drug generalisation. |
| `TranSiGen/` | [TranSiGen](https://github.com/myzhengSIMM/TranSiGen) `e7582e77` | MIT | Stage 2 single-dose entrant | Transcriptome generative model; random initialization and Morgan2048 under ADR 0074; Stage 1 provenance retained. |
| `XPert/` | author-provided (version unknown) | MIT | reference | Stage 1 XPert library; `local/` holds the routing index, logs, metadata, metrics and figure-reproduction code. |

Nine items historically completed Stage 2 setup; their run cards are retained in
Git history. Five were historically **reference-only** (`WAVE`, `ExPO`, `State`, `TranSiGen`,
`XPert`). ADR 0074 now admits TranSiGen to the single-dose 5k comparison,
without using retained pretrained assets or changing Stage 3. The current successor comparison is [04_model_view](../stage2/04_model_view/README.md),
using the released [HVG-DEG expression views](../stage2/03_view/README.md) inputs. `pert_data/`,
`drug_emb_model/`, `disease_data/` are the sibling external categories under the
same contract; `pert_model/` is the model-code one, so its `source/` trees are
the upstream repositories themselves plus whatever weights/data the authors
publish separately.

## Addition on 2026-11-05

`PerturbDiff/` contains the MIT-licensed
[DeepGraphLearning/PerturbDiff](https://github.com/DeepGraphLearning/PerturbDiff)
source at commit `31be89cb8cce5dddef0782bac29d5c748ce55931`,117 immutable
files/1,203,863bytes. Its `author_log.json`, source manifest and offline verifier
record provenance; no upstream checkpoint or dataset was downloaded.
Historical native setup proved 2,000-gene synthetic native training, save/reload
and DDIM sampling in an isolated overlay. Its setup and tenv adapter/runtime
records are retained in Git history; native setup did not establish formal
benchmark admission. The fourteen-item inventory above describes that historical
setup status.

## Excluded preprints (not part of this repository)

The old `download/manual.md` exclusion table, kept here because it is the
standing admission record:

| model | reason |
|---|---|
| X-Pert, DORA, BaiZe, PerturbLDM | preprints; excluded from Stage 2 setup and removed from `pert_model/` |
| WAVE | reference-only; no dose input |
| ExPO | reference-only; 978-only published scope |
| DeepICER, MiTCP | 978-only and no LICENSE |
| CellPert, scNODE, trVAE | no dose input |
| PertFlow | no LICENSE / no README |

The full revised-gate survey and model admission evidence are retained in Git
history with the removed Stage 2 audit and setup trees.

## Manual downloads — after downloading

Follow each item's `script/download.md`, retain downloaded bytes under `source/`,
and run its `script/verify.sh` for offline source verification. The historical
model setup checks are available in Git history; this checkout has no current
model setup or benchmark runner.

Items that could not be fetched automatically (FigShare returns HTTP 403 / a WAF
challenge from this WSL host; Google Drive is browser-only) carry the browser
instructions in their own `script/download.md`: WAVE, Squidiff, scVIDR, DECANT
(no public data), plus the author-provided XPert library. PerturbNet's assets
were fetched through the `hf-mirror.com` mirror of
`cyclopeta/PerturbNet_reproduce`; CMonge's demo data is bundled in its upstream
repository; PRnet, CPA, chemCPA and biolord carry their native assets in the
positions recorded in their `author_log.json` `sources[]`.
