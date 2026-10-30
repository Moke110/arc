# Download — chemCPA

**method: repo + fetched asset** — upstream clone (MIT) plus the authors' public
model archive on Backblaze B2.

## Repository

```bash
git clone https://github.com/theislab/chemCPA pert_model/chemCPA/source
git -C pert_model/chemCPA/source checkout 43e830eb0958c54e4aa64442c17ec0fed19b3f15
# vendored read-only: the nested .git directory was removed on 2026-10-18
```

## Weights and data

The released weights are a public archive, not a browser-only download:

```bash
wget https://f003.backblazeb2.com/file/chemCPA-models/chemCPA_models.zip \
  -O pert_model/chemCPA/source/source/chemCPA_models.zip
cd pert_model/chemCPA
unzip source/source/chemCPA_models.zip -d source/weights/
```

`source/weights/chemCPA_models/` holds 12 raw ComPert state dicts (L1000
977/2,000 genes, 188 drugs). `source/data/combo_sciplex_prep_hvg_filtered.h5ad`
is the 5,000-HVG ComPert dataset kept for future retraining; the Stage 2 smoke
is load-only and does not need it.

## Cleanup applied

`weights/__MACOSX/` (macOS archive metadata) and the dangling
`project_folder` symlink to `/nfs/homedirs/hetzell/hdd/project_chemCPA` were
removed in ticket 04; they were not download content. SHA256s are in
`author_log.json`.
