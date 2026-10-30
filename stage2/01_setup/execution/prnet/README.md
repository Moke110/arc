# PRnet setup

**Status:** setup-only (978-gene LINCS demo; not a Stage 2 evaluation result)

- source: `pert_model/PRnet/source/` (commit `f19174bd`), Apache-2.0
- weights: `pert_model/PRnet/source/checkpoint/lincs_best_epoch_all.pt` (0.55 M params, 978 genes)
- data: `pert_model/PRnet/source/dataset/demo.h5ad` (681 × 978, `demo_split`)
- env: `arc-stage2-bringup-prnet` (Python 3.9, CPU torch)
- command:

```bash
bash stage2/01_setup/execution/prnet/run.sh
```

Expected output: `output/run.log`. The 2026-10-19 run wrote
`demo_split_{y_pre,x_true,y_true,cov_drug}_array.csv` (the 661 × 978 predictions,
48 MB); deleted in ticket 06, a rerun recreates them.

**Result (2026-10-19 re-run):** passed. The test split produced **661 × 978**
predictions; the checkpoint was loaded on CPU via a small `torch.load`
`map_location` wrapper (`run_prnet.py`) because upstream calls `torch.load`
without it. No upstream code was edited.

Caveat: the 978-landmark matrix is allowed **only** here. It must not be read by
any Stage 2 evaluation or compared with whole-transcriptome models.

**5,000-gene path (not exercised).** The clone also ships
`checkpoint/sciplex_best_epoch_all.pt` and `test_sciplex.py`, but the latter reads
`dataset/Sci_Plex.h5ad`, which upstream does not distribute. Running the 5,000
path therefore needs that processed Sci-Plex matrix first.
