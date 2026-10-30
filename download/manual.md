# Manual downloads — Stage 2 setup

Only files that **could not** be fetched automatically are listed here. FigShare
returns HTTP 403 / a WAF challenge from this WSL host, so those need a browser
download (or a different network).

Already downloaded automatically:

- `pert_model/perturbnet/assets/` — 2.1 GB of LINCS/sci-Plex checkpoints +
  examples (via `hf-mirror.com` mirror of `cyclopeta/PerturbNet_reproduce`)
- `pert_model/cmonge/data/` — bundled demo data
- the nine deployed model repositories under `pert_model/`
- `pert_model/prnet`, `cpa`, `chemcpa`, `biolord` — Stage 2 native assets

Stage 2 setup excludes preprints (X-Pert, DORA, BaiZe, PerturbLDM); their code and
downloads are no longer part of this repository.

---

## 1. WAVE — pretrained model + preprocessed datasets (reference-only)

WAVE is not deployed under the revised gates (no dose input), but its code is kept
for comparison. The repository ships no weights; the paper's pretrained model and
preprocessed LINCS/sci-Plex h5ad files are on FigShare.

**Source:** <https://figshare.com/articles/dataset/Pretrained_WAVE_Model_and_Datasets_for_Drug-Induced_Gene_Expression_Prediction_and_TCGA-Based_Drug_Screening/29650832>
(DOI `10.1038/s41467-026-77645-3`, FigShare `29650832`)

1. Open the page in a browser and use **Download all**.
2. Put the archive(s) in `pert_model/wave/source/`.
3. Extract the checkpoint to `pert_model/wave/weights/`; the training code
   expects a `best_model.pth`-style state dict for `wave.model.WAVE`.
4. Extract the preprocessed `.h5ad` files to `pert_model/wave/data/`; each needs
   `X`, `layers['unpert_expr']`, and `obs['smiles']` (what
   `wave/load_dataset.py` reads).
5. Record SHA256s in `pert_model/wave/VENDORED.txt`.

Fallback raw data: LINCS Level 3 (<https://clue.io/releases/data-dashboard>),
sci-Plex3 GEO `GSM4150378`.

---

## 2. Squidiff — released checkpoints (deployed)

The repository ships code but no weights; the paper's checkpoints are on FigShare.

**Source:** <https://doi.org/10.6084/m9.figshare.27948633> (FigShare `27948633`)

1. Open the DOI in a browser and download all checkpoint files.
2. Put them in `pert_model/squidiff/weights/`.
3. Keep the original names; `sample_squidiff.sampler(model_path=...)` loads one
   `.pt` checkpoint.
4. Record SHA256 and file names in `pert_model/squidiff/VENDORED.txt`.

Companion reproducibility repo:
<https://github.com/siyuh/Squidiff_reproducibility>.

---

## 3. scVIDR — `data.zip` (deployed; includes pretrained models)

The repository has code only; the figure/notebook data and pretrained
single-/multi-dose models are in the author's Google Drive archive.

**Source:** <https://drive.google.com/file/d/11fzDbp0B19Dy47MtD742Jl4Hz2bdSiiq/view?usp=sharing>

1. Download `data.zip`.
2. Put it in `pert_model/scvidr/source/` and extract so the notebooks see
   `pert_model/scvidr/data/`.
3. The archive contains `nault2021_singleDose.h5ad`,
   `nault2021_multiDose.h5ad`, and the reported `.pt` model directories.
4. Record the archive SHA256 in `pert_model/scvidr/VENDORED.txt`.

The 2,000-gene synthetic smoke in `stage2/01_setup/scvidr/` does **not** need
this archive; it is only required for the published liver dose-response results.

---

## 4. DECANT — no public data file

DECANT (deployed) explicitly excludes processed cell bags, teacher prototypes,
and checkpoints from its repository. There is no public download. The setup smoke
uses a synthetic 2,000-gene forward; reconstructing matched treated/control cell
bags is a separate data-engineering task.

---

## 5. Not part of this repository

| model | reason |
|---|---|
| X-Pert, DORA, BaiZe, PerturbLDM | preprints; excluded from Stage 2 setup and removed from `pert_model/` |
| WAVE | reference-only; no dose input |
| ExPO | reference-only; 978-only published scope |
| DeepICER, MiTCP | 978-only and no LICENSE |
| CellPert, scNODE, trVAE | no dose input |
| PertFlow | no LICENSE / no README |

---

## 6. After downloading

Run the checks in `stage2/01_setup/`. For WAVE the minimum real-data smoke is one
fold of the preprocessed LINCS data with a single epoch; for Squidiff, load one
released checkpoint and run `sample_squidiff.sampler` on a small test h5ad; for
scVIDR, train the multi-dose model on the `data.zip` AnnData and compare
`predict(continuous=True, doses=[...])`. None of these needs a code change inside
`pert_model/<model>/`.
