# PerturbNet setup smoke

**Status:** checkpoint/architecture smoke (no TF1.15 cell-VAE inference)

- source: `pert_model/PerturbNet/source/` (commit `e1375b2c`), GPL-3.0
- assets: `pert_model/PerturbNet/source/assets/` — LINCS 978 and sci-Plex 5,087 cINN
  checkpoints, ChemicalVAE, cell VAEs; see `pert_model/PerturbNet/source/ASSETS.md`
- env: `arc-stage2-bringup-perturbnet` (Python 3.10, CPU torch, numpy/pandas/rdkit)
- command:

```bash
bash stage2/01_setup/execution/perturbnet/run.sh
```

Expected output: `output/meta.json`, `output/run.log`. The two
`pred_latent_*.npy` arrays were deleted in ticket 06; a rerun recreates them.

Caveat: the stable PerturbNet package needs Python 3.7 + TensorFlow 1.15 for
end-to-end cell-VAE inference. This smoke loads the released cINN + ChemicalVAE
weights into the in-repo architecture and runs a latent forward; it proves
checkpoint compatibility, not a CM prediction.

**Result (2026-10-19):** passed. LINCS and sci-Plex cINN checkpoints loaded with all 816 tensors matched; ChemicalVAE loaded (4,030,982 params); latent forward produced (2, 10) for both checkpoints.
