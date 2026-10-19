# AGENTS.md — stage2/01_setup

One directory per deployed model. Each run card is a README plus a `run.sh`;
outputs go to `<model>/output/` (gitignored).

- `prnet/` — LINCS 978 demo (`test_demo.py`); **setup-only**, not a Stage 2
  evaluation artifact.
- `cpa/` — Kang checkpoint + `kang_normalized_hvg.h5ad`.
- `chemcpa/` — load-only forward with the public L1000 977/2,000-gene checkpoint.
- `biolord/` — small from-scratch training on `GSM_new.h5ad`.
- `perturbnet/` — released LINCS/sci-Plex cINN + ChemicalVAE load and latent forward.
- `scvidr/` — synthetic 2,000-gene train and continuous-dose prediction.
- `squidiff/` — 2,000-gene drug-structure forward.
- `cmonge/` — official sciplex demo train + evaluate.
- `decant/` — synthetic 2,000-gene forward with two dose values.

Run from the repository root:

```bash
bash stage2/01_setup/<model>/run.sh
```

Every run writes `output/run.log` and at least one artifact. Do not edit
`pert_model/<model>/` to make a run work; put adapters under the model's setup
directory. Preprint models (X-Pert, DORA, BaiZe, PerturbLDM) are not part of this
stage.
