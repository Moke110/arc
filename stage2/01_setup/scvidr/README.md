# scVIDR setup smoke

**Status:** synthetic 2,000-gene train + continuous-dose prediction

- source: `pert_model/scvidr/` (commit `bcb6a96c`), MIT
- weights: none published; the smoke trains a small model from scratch
- data: synthetic AnnData (3 cell types × doses 0/1/10/100 × 20 cells × 2,000 genes)
- env: `arc-stage2-bringup-scvidr` (Python 3.8, scvi-tools 0.13.0,
  pytorch-lightning 1.3.8, torch 2.4.1 CPU, scanpy 1.9.8, geomloss 0.2.5)
- command:

```bash
bash stage2/01_setup/scvidr/run.sh
```

Expected output: `output/pred_doses.npz`, `output/meta.json`, `output/run.log`.

Caveat: the continuous path is a log-linear latent interpolation
(`delta * log1p(d) / log1p(max(doses))`), not a learned dose-response curve; the
shipped CLI rejects doses absent from `obs`.

**Result (2026-10-19):** passed. Trained on 180 synthetic cells × 2,000 genes for 2 epochs; continuous predictions returned (15, 2,000) for doses 1, 10 and 100 µM.
