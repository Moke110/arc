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

`manifest.sha256` (per-file hashes over `source/` only) exists for every item — 442 GB hashed
across the four categories, 27/27 verified. `script/verify.sh` checks it offline, and a
missing manifest is a failure. Run `script/verify_all.sh` from the repository root to check
every item.

## The fourteen 项

| 项 | upstream | licence | stage | what it is |
|---|---|---|---|---|
| `PRnet/` | [PRnet](https://github.com/Perturbation-Response-Prediction/PRnet) `f19174bd` | Apache-2.0 | deployed | PRnet reaction/response predictor; ships its own 978 LINCS + 5,000 Sci-Plex checkpoints. |
| `CPA/` | [CPA](https://github.com/theislab/CPA) `fbd7c025` | BSD-3-Clause | deployed | Compositional Perturbation Autoencoder; Kang + Norman2019 pretrained models on the authors' Drive. |
| `chemCPA/` | [chemCPA](https://github.com/theislab/chemCPA) `43e830eb` | MIT | deployed | Chemistry-aware CPA; public 977/2,000-gene L1000 ComPert weights, not retrained on CM data. |
| `biolord/` | [biolord](https://github.com/nitzanlab/biolord) `b7688790` | BSD-3-Clause | deployed | Disentangled representation model; no public weights (per-dataset training). |
| `PerturbNet/` | [PerturbNet](https://github.com/welch-lab/PerturbNet) `e1375b2c` | GPL-3.0 | deployed | Conditional INN + ChemicalVAE; 5,087-gene sci-Plex checkpoint with a raw continuous dose covariate. |
| `scVIDR/` | [scVIDR](https://github.com/BhattacharyaLab/scVIDR) `bcb6a96c` | MIT | deployed | VAE with single-/multi-dose modes; 5,000 HVG, log-linear latent dose interpolation. |
| `Squidiff/` | [Squidiff](https://github.com/siyuh/Squidiff) `abdfc27d` | BSD-3-Clause | deployed | Diffusion model; configurable gene space, drug input is rFCFP × log10(dose+1). |
| `CMonge/` | [Conditional-Monge](https://github.com/AI4SCR/Conditional-Monge) `195f5cd1` | MIT | deployed | Conditional Monge-gap transport; RDKit structure + categorical `drug-dose`, bundled demo data. |
| `DECANT/` | [DECANT](https://github.com/bliulab/DECANT) `e50642a9` | MIT | deployed | Mechanism/context-decoupled drug representation; dose is context to suppress, no shipped data/weights. |
| `WAVE/` | [WAVE](https://github.com/ZJUFanLab/WAVE) `e27e453b` | MIT | reference | GeneDrugFusion; 978 default, **no dose input** → fails revised G2. |
| `ExPO/` | [ExPO](https://github.com/MLBC-lab/ExPO) `7e88cbae` | MIT | reference | Exposure-aware model; **978 only** in the published scope → fails revised G1. |
| `State/` | [state](https://github.com/ArcInstitute/state) `9bbfe78a` | CC BY-NC-SA 4.0 / State Model Non-Commercial | reference | Arc State transition model (SE-600M + ST-HVG-Tahoe); Tahoe one-hot limits cold-drug generalisation. |
| `TranSiGen/` | [TranSiGen](https://github.com/myzhengSIMM/TranSiGen) `e7582e77` | MIT | reference | Transcriptome generative model; kept for Stage 1 provenance. |
| `XPert/` | author-provided (version unknown) | MIT | reference | Stage 1 XPert library; `local/` holds the routing index, logs, metadata, metrics and figure-reproduction code. |

Nine items are **deployed** in Stage 2 setup (run cards under
`../stage2/01_setup/`); five are **reference-only** (`WAVE`, `ExPO`, `State`,
`TranSiGen`, `XPert`) and have no Stage 2 run card. `pert_data/`,
`drug_emb_model/`, `disease_data/` are the sibling external categories under the
same contract; `pert_model/` is the model-code one, so its `source/` trees are
the upstream repositories themselves plus whatever weights/data the authors
publish separately.

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

The full revised-gate survey is the relocated Stage 2 audit record
[`../stage2/audit/model_survey_2026-10-19.md`](../stage2/audit/model_survey_2026-10-19.md).

## Manual downloads — after downloading

The deleted combined note's generic procedure, kept here:

Run the checks in `../stage2/01_setup/`. For WAVE the minimum real-data smoke is
one fold of the preprocessed LINCS data with a single epoch; for Squidiff, load
one released checkpoint and run `sample_squidiff.sampler` on a small test h5ad;
for scVIDR, train the multi-dose model on the `data.zip` AnnData and compare
`predict(continuous=True, doses=[...])`. None of these needs a code change inside
an item's `source/` tree.

Items that could not be fetched automatically (FigShare returns HTTP 403 / a WAF
challenge from this WSL host; Google Drive is browser-only) carry the browser
instructions in their own `script/download.md`: WAVE, Squidiff, scVIDR, DECANT
(no public data), plus the author-provided XPert library. PerturbNet's assets
were fetched through the `hf-mirror.com` mirror of
`cyclopeta/PerturbNet_reproduce`; CMonge's demo data is bundled in its upstream
repository; PRnet, CPA, chemCPA and biolord carry their native assets in the
positions recorded in their `author_log.json` `sources[]`.
