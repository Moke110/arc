# Download — WAVE

**method: manual** — reference-only. The repository ships no weights; the paper's
pretrained model and preprocessed LINCS/sci-Plex h5ad files are on FigShare.

## Repository (automatic)

```bash
git clone https://github.com/ZJUFanLab/WAVE pert_model/WAVE/source
git -C pert_model/WAVE/source checkout e27e453beb48d0ff912585060b5a09cb04335262
# vendored read-only: the nested .git directory was removed on 2026-10-19
```

## 1. WAVE — pretrained model + preprocessed datasets (reference-only)

WAVE is not deployed under the revised gates (no dose input), but its code is kept
for comparison. The repository ships no weights; the paper's pretrained model and
preprocessed LINCS/sci-Plex h5ad files are on FigShare.

**Source:** <https://figshare.com/articles/dataset/Pretrained_WAVE_Model_and_Datasets_for_Drug-Induced_Gene_Expression_Prediction_and_TCGA-Based_Drug_Screening/29650832>
(DOI `10.1038/s41467-026-77645-3`, FigShare `29650832`)

1. Open the page in a browser and use **Download all**.
2. Put the archive(s) in `pert_model/WAVE/source/`.
3. Extract the checkpoint to `source/weights/`; the training code
   expects a `best_model.pth`-style state dict for `wave.model.WAVE`.
4. Extract the preprocessed `.h5ad` files to `source/data/`; each needs
   `X`, `layers['unpert_expr']`, and `obs['smiles']` (what
   `wave/load_dataset.py` reads).
5. Record SHA256s in `author_log.json` sources[] (`VENDORED.txt` is retired).

Fallback raw data: LINCS Level 3 (<https://clue.io/releases/data-dashboard>),
sci-Plex3 GEO `GSM4150378`. The FigShare assets were not downloaded on this host;
the synthetic 978-gene smoke in the deleted VALIDATION.md used random weights.
