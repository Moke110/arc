# AGENTS.md — pert_model
Drug-perturbation model code and weights as external 项 (ADR 0030):
each has tracked `script/` + `author_log.json`, and gitignored `source/`
(upstream clone + downloaded assets) and `local/` (declared local additions).
Historically set up (run cards are retained in Git history):

- `PRnet/` — PRnet (Apache-2.0); 978 LINCS + 5,000 Sci-Plex checkpoints.
- `CPA/` — CPA (BSD-3-Clause); Kang + Norman2019 checkpoints, 5,000 genes.
- `chemCPA/` — chemCPA (MIT); 977/2,000-gene L1000 ComPert weights, not retrained.
- `biolord/` — biolord (BSD-3-Clause); no public weights, per-dataset training.
- `PerturbNet/` — PerturbNet (GPL-3.0); 5,087-gene sci-Plex checkpoint + dose.
- `scVIDR/` — scVIDR (MIT); 5,000 HVG, log-linear latent dose interpolation.
- `Squidiff/` — Squidiff (BSD-3-Clause); configurable genes, rFCFP × dose scaling.
- `CMonge/` — Conditional Monge Gap (MIT); RDKit structure + categorical `drug-dose`, bundled demo data.
- `DECANT/` — DECANT (MIT); configurable genes, dose context, no data/weights.
- `PerturbDiff/` — native setup and tenv admission evidence are retained in Git history.
- `TranSiGen/` — TranSiGen (MIT); random-initialized native Stage 2 single-dose 5k entrant under ADR 0074; retained pretrained assets are not used.

Reference-only:

- `WAVE/` — WAVE (MIT); 978 default, **no dose input** → fails revised G2.
- `ExPO/` — ExPO (MIT); **978 only** in published scope → fails revised G1.
- `State/` — Arc State (non-commercial weights); Tahoe one-hot limits cold drug.
- `XPert/` — author-provided Stage 1 XPert library; no public download.

Layout: ADR 0030. Category overview, exclusion table and the manual-download
procedure: `README.md`. Every item carries `manifest.sha256` (453.0 GB hashed,
27/27 verified); `script/verify.sh` checks it offline.

Current Stage 2 comparison adapters are routed from [the Stage 2 index](../stage2_CarPS-VAE/AGENTS.md). Verify retained external sources with each item’s `script/verify.sh`.
