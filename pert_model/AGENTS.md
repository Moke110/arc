# AGENTS.md — pert_model

Drug-perturbation model code and weights. Upstream clones are vendored read-only;
their `VENDORED.txt` records commit, licence, and assets. Large data/weights are
gitignored.

## Deployed models (Stage 2 setup)

- `prnet/` — PRnet (Apache-2.0); 978 LINCS + 5,000 Sci-Plex checkpoints.
- `cpa/` — CPA (BSD-3); Kang checkpoint, 5,000 genes.
- `chemcpa/` — chemCPA (MIT); public weights are 977/2,000-gene L1000 ComPert models.
- `biolord/` — biolord (BSD-3); no public weights, per-dataset training.
- `perturbnet/` — PerturbNet (GPL-3.0); 5,087-gene sci-Plex checkpoint + dose.
- `scvidr/` — scVIDR (MIT); 5,000 HVG, continuous latent dose interpolation.
- `squidiff/` — Squidiff (BSD-3); configurable gene space, rFCFP × dose scaling.
- `cmonge/` — Conditional Monge Gap (MIT); RDKit structure + categorical
  `drug-dose` conditions, per-dataset HVG.
- `decant/` — DECANT (MIT); configurable genes, dose context, no data/weights.

Run cards and results: `../stage2/01_setup/`.

## Reference-only (not deployed)

- `wave/` — WAVE (MIT); 978 default, **no dose input** → fail revised G2.
- `expo/` — ExPO (MIT); **978 only** in published scope → fail revised G1.
- `state/` — Arc State (non-commercial weights); Tahoe one-hot limits cold drug.
- `TranSiGen/` — TranSiGen; kept for Stage 1 provenance.
- `xpert/` — Stage 1 XPert library.

## Excluded preprints

DORA, BaiZe, and PerturbLDM are preprints with no Stage 2 vendor clone or run
card. X-Pert is also excluded from Stage 2, but the Stage 1 XPert library
remains at `xpert/` and is not a Stage 2 entrant. All four remain documented in
`SURVEY_2026-10-19_v2.md` and `BIOMNI_SURVEY_2026-10-19.md`.

## Documentation

- `SURVEY_2026-10-19_v2.md` — revised-gate survey (current admission record).
- `SURVEY_2026-10-19.md` — initial four-gate survey (historical).
- `VALIDATION.md` — earlier 2,000-gene smoke log; Stage 2 setup has its own run
  cards under `../stage2/01_setup/`.
- `../download/manual.md` — manual downloads (WAVE, Squidiff, scVIDR).

## Data manifest

`STAGE2_DATA.md` records the h5ad/zip provenance and SHA256 for the Stage 2
setup assets.
