# Squidiff setup smoke

**Status:** 2,000-gene drug-structure forward with random weights and a real
rFCFP/dose encoding

- source: `pert_model/squidiff/` (commit `abdfc27d`), BSD-3
- weights: none released; FigShare checkpoints are a manual download, so the
  network is randomly initialised in this smoke
- data: synthetic 2,000-gene expression batch; the drug condition is computed
  from a real SMILES string (dexamethasone) at doses 1 and 100 µM
- env: `arc-stage2-bringup-squidiff` (Python 3.13, CPU torch, numpy/rdkit)
- command:

```bash
bash stage2/01_setup/squidiff/run.sh
```

Expected output: `output/pred.npy`, `output/meta.json`, `output/run.log`.

Caveat: dose enters through an rFCFP fingerprint multiplied by
`log10(dose+1)`, not as a separate learned dose curve. The drug encoding is
real; the network weights are random, so the output values are not biological
predictions.

**Result (2026-10-19 re-run):** passed. 52,559,884-parameter model; the forward
pass used `Squidiff.Drug_dose_encoder`'s encoding (rFCFP × `log10(dose+1)`) from
dexamethasone SMILES at 1 and 100 µM, and produced (2, 2,000).
