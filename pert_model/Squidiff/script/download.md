# Download — Squidiff

**method: manual** — deployed model. The repository ships code but no weights; the
paper's checkpoints are on FigShare.

## Repository (automatic)

```bash
git clone https://github.com/siyuh/Squidiff pert_model/Squidiff/source
git -C pert_model/Squidiff/source checkout abdfc27d84947dcccd745d1067c0840a41d32eb8
# vendored read-only: the nested .git directory was removed on 2026-10-19
```

`LICENSE` still contains unresolved merge-conflict markers; both sides are
BSD-3-Clause (azizilab 2022, Siyu He 2024).

## 2. Squidiff — released checkpoints (deployed)

The repository ships code but no weights; the paper's checkpoints are on FigShare.
FigShare returns HTTP 403 / a WAF challenge from this WSL host, so this is a
browser download.

**Source:** <https://doi.org/10.6084/m9.figshare.27948633> (FigShare `27948633`)

1. Open the DOI in a browser and download all checkpoint files.
2. Put them in `source/weights/`.
3. Keep the original names; `sample_squidiff.sampler(model_path=...)` loads one
   `.pt` checkpoint.
4. Record SHA256 and file names in `author_log.json` sources[] (`VENDORED.txt`
   is retired).

Companion reproducibility repo:
<https://github.com/siyuh/Squidiff_reproducibility>. The checkpoints were not
downloaded on this host; the synthetic 500-gene forward pass in the deleted
VALIDATION.md used randomly initialised weights.
