# Download — CPA

**method: repo + manual asset** — upstream clone (BSD-3-Clause) plus the authors'
public Google Drive "Datasets and Pre-trained models" folder.

## Repository

```bash
git clone https://github.com/theislab/CPA pert_model/CPA/source
git -C pert_model/CPA/source checkout fbd7c0250edc23eff003a10c99655579c53afd63
# vendored read-only: the nested .git directory was removed on 2026-10-18
```

The clone's code, docs, examples and tests stay at their clone-relative position
under `source/`.

## Datasets and weights (manual browser download)

The Drive folder <https://drive.google.com/drive/folders/1yFB0gBr72_KLLp1asojxTgTqgz6cwpju>
does not serve an unattended download, so this part is **manual**:

1. Open the folder in a browser and download
   `Pretrained_Models-20261018T124453Z-1-001.zip`; put it under
   `pert_model/CPA/source/source/` (upstream-style `source/` subdirectory of the clone)
   and record its SHA256 (see the author log).
2. Download the two dataset h5ads and put them under `pert_model/CPA/source/data/`:
   `kang_normalized_hvg.h5ad`, `Norman2019_normalized_hvg.h5ad`.
3. Extract the zip so the models land under `pert_model/CPA/source/weights/Pretrained_Models/`:

```bash
cd pert_model/CPA
unzip source/source/Pretrained_Models-20261018T124453Z-1-001.zip -d source/weights/
```

## Notes

Deployed Stage 2 model; run card `stage2/01_setup/execution/cpa/` uses the Kang checkpoint
(5,000 genes). SHA256s for both h5ads and the zip are recorded in
`author_log.json`. `manifest.sha256` (ticket 09) covers `source/` only.
