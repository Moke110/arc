# Download — PRnet

**method: repo** — upstream clone, pinned commit, Apache-2.0. No separately
downloaded asset: the released checkpoints ship inside the repository.

## Repository

The upstream clone that produced `source/`:

```bash
git clone https://github.com/Perturbation-Response-Prediction/PRnet \
  pert_model/PRnet/source
git -C pert_model/PRnet/source checkout f19174bde3ed2633f54c7831799cc38c4ffc7a0d
# vendored read-only: the nested .git directory was removed on 2026-10-18
```

Everything the upstream repository ships stays at its clone-relative position
under `source/`; nothing was edited.

## Assets

`checkpoint/` (inside the clone) holds the released state dicts verified in the
Stage 2 setup smoke:

| file | bytes | shape |
|---|---:|---|
| `checkpoint/lincs_best_epoch_all.pt` | 2,210,427 | 978-landmark LINCS |
| `checkpoint/sciplex_best_epoch_all.pt` | 8,420,411 | 5,000-HVG Sci-Plex |

`dataset/`, `data/`, `models/`, `trainer/`, `preprocessing/`, `figure/`, `img/`
are repository content and were not modified.

## Notes

Deployed Stage 2 model; run card `stage2/01_setup/execution/prnet/`. The per-file hashes
arrive with `manifest.sha256` in ticket 09; `script/verify.sh` is a declared
no-op until then.
