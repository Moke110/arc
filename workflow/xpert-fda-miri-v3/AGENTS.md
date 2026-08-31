# AGENTS.md

The **all-CM finetune MIRI screening workflow** (`workflow/xpert-fda-miri-v3/`;
see docs/adr/0017): the 11-step sibling of the finalized
`workflow/xpert-fda-miri-v1/` (ADR 0007) whose **only substantive change is the
XPert finetune data**. Instead of v1's frozen DToxS-only `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`,
steps 05/06/07 train on the workflow-built **union h5ad of all four usable CM
datasets** (`dtoxs_primary` from the LINCS LDS-1593 well-level package,
`dtoxs_ipsc` from GSE217423, `gse244740` ScreenSeq, `gse262419` TempO-Seq,
step 03; ADR 0017), under the **per-plate normalization** of ADR 0016
(matched control + QN per `(dataset, plate)`). gse247541 is excluded (6-day
exposure outside the author time map, 2 in-HG drugs, 12 treated samples).
Since ADR 0015 the DToxS arms carry **real-concentration dose bins** (not v1's
single-dose convention) with a single-cell encoding, and the benchmark test
sets are **byte-identical to v2's** (fold columns imported verbatim for the
primary arm; the two concentration-response datasets are train-only). v2
stays authoritative and untouched. There is **no agentic report stage** (the
former `12_report` was removed per spec 1; v1 archives its report as a
non-numbered `report/` stage).

## Scope and provenance (author + CM-dataset data only)

- Screens **only** the FDA drugs already present in the L1000 embedding table
  (exact canonical-SMILES match, `drug_idx < 8981`; 1,004 of 3,127 drugs).
  Because every screened drug uses the author-reference HG embedding directly,
  **no HG training / low-LR finetune / FDA-extension inference is ever run**.
- **Drug info** comes only from the workflow-owned `input/` files (author L1000
  SMILES/HG tables, precomputed UniMol array, inherited FDA query list) and
  the per-dataset files (incl. the DToxS small-molecule metadata for SMILES).
  The workflow never queries an external drug database and never regenerates
  drug features from SMILES.
- **Gene info** comes only from the XPert-author
  `pert_model/xpert/metadata/l1000_gene_info_978.csv` and each dataset's own
  gene index. No external gene databases.
- `pert_data/` is **read-only**: step 03 never modifies a source file; the
  union h5ad + maps + manifest are regenerated workflow artifacts.
- The **"no external database" provenance rule is amended for the dose field
  only** (ADR 0015): the primary arm's finetune concentrations are cited
  (paper Table 1 KIs / GSE217423 series-matrix shared non-KIs / literature
  clinical Cmax for the 8 primary-only drugs) in
  `03_cm_finetune_data/primary_concentrations.py`.

## Structure (11-step layout)

- `01_drug_embed/` — merged drug list (`drug_list.py`, 1,004 drugs) + author-reference HG embeddings (`drug_embeddings.py`, byte-parity gate). Identical to v1.
- `02_cm_baseline/` — consensus baseline from the DToxS-primary control wells of the released LINCS LDS-1593 package (`pert_data/dtoxs_primary/Data/Level3`); byte-identical to v1's baseline.
- `03_cm_finetune_data/` — **the v3 change**: builds one finetune h5ad per dataset (`output/dtoxs_primary.h5ad` 662 rows, `output/dtoxs_ipsc.h5ad` 296 rows, `output/gse244740.h5ad` ~2,424 rows, `output/gse262419.h5ad` ~936 rows) and concatenates them into the union `output/cm_union.h5ad` (~4,318 rows), per-plate normalized (ADR 0016), plus per-dataset gene maps (`output/gene_maps/`), the compound→drug-index manifest (`output/drug_manifest.csv`), a schema gate vs v1's `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`, the delta report vs v2's union (`output/finetune_delta.md`), and the composition report `output/finetune_report.html`.
- `04_xpert_pretrain/` — the workflow-owned pretrain unit (weight + `config/l1000.yaml`), carried over untouched.
- `05_xpert_finetune_bench/` — union benchmark on the **imported v1.1 folds** (ADR 0014/0015/0017): `create_cm_splits.py` imports v1.1's `warm-start_fold-1..10` + `cold-drug_fold-1..10` verbatim for primary rows (matched by well id, byte-identical test sets), assigns iPSC rows under the ADR 0015 isolation rules, and assigns **gse244740/gse262419 rows train-only** (never test; cold-drug fold test drugs excluded from all four datasets' train) + `cm_full` (all-train), with four isolation gates; `xpert_finetune_bench.py` benches the 20 folds at **4 epochs** and — as the **first column** (`pretrain`, amber) — zero-shot-evaluates the step-04 L1000 pretrain weight on all 20 fold test sets (no finetuning) under the workflow-owned `config/cm_pretrain_eval.yaml`. Comparable column-by-column with v2 (identical test sets).
- `06_xpert_finetune/` — all-train union finetune (`xpert_finetune.py`, `config/cm.yaml`) -> the final prediction model `output/cm_cm_full/`.
- `07_xpert_predict/` — perturbed-transcriptome prediction + z-score (`perturb_transcriptome.py`, `zscore.py`); predicts at dose_idx 0 / time_idx 2 (the v1 exposure, unchanged by the real-dose finetune).
- `08_miri_phenotype/` — integrated MIRI reference gene sets + rank reference.
- `09_miri_efficacy_scores/` — the three efficacy scorers + z-components + plots.
- `10_filter/` — WCS<0 + >=6/8 critical-pathway filters.
- `11_integrate/` — integrated score over all drugs, ranked candidates, deterministic report.
- `input/` — frozen assets (author HG reference, UniMol array, idx→SMILES, FDA drug list, curated `MIRI_pathways.json`), hard-linked from v1.
- `_lib/utils.py` — workflow-root-aware paths; `CM_DATA_ROOT`, `CM_DATASETS` (all five, for the data-lake audit), `FINETUNE_DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse262419"]`, `N_L1000_DRUGS = 8981`, `AUTHOR_TIME_MAP = {3:0, 6:1, 24:2}`, `L1000_PRETRAIN_H5AD`, `AUTHOR_HG_REFERENCE`, `V11_DTOXS_LINCS_H5AD` (fold source), `V2_CM_UNION_H5AD` (delta source).
- `_lib/gene_mapping.py` — **per-dataset** CM→L1000 landmark maps (symbols / Ensembl / SYMBOL_entrezid), legacy-symbol aware, `corrected_gene_name` never a match key; gse262419's TempO-Seq platform-missing landmarks {WRB, FAM57A} are zeroed.
- `_lib/drug_resolution.py` — compound→L1000 drug index (canonical-SMILES match, else author `cmap_name` name match) + the auditable manifest.
- `_lib/verify_genes.py` — the extended per-dataset gene gates (G1 completeness, G2 alias verification; G3 union-column audit for the finetune datasets only), writes `_lib/gene_mapping_report.md`.

## The union finetune data (step 03, ADR 0017)

- The finetune scope is **all four usable CM datasets** (`dtoxs_primary` from
  the LINCS LDS-1593 well-level package, `dtoxs_ipsc` from GSE217423,
  `gse244740` ScreenSeq, `gse262419` TempO-Seq); gse247541 is excluded
  (6-day exposure, 2 in-HG drugs, 12 treated samples) but stays in the data
  lake + drug-overlap audit. The union is **~4,318 rows / ~180 unique
  drug_idx**: dtoxs_primary 662, dtoxs_ipsc 296, gse244740 ~2,424 (38 drugs),
  gse262419 ~936 (101 in-HG drugs of 464 chemicals — the TempO-Seq panel is
  mostly environmental/toxicology compounds absent from the author table, so
  gse244740, not gse262419, dominates the union).
- Each dataset is preprocessed independently: 978 landmarks → CPM → log2(CPM+1)
  → per-`(dataset_id, plate)` quantile alignment to the L1000 pretrain pooled
  distribution (one monotonic transform per plate, before matched-control
  computation; QN reference cached at `output/qn_reference.npy`) → **per-plate
  matched control** (each sample's baseline = the mean of the control wells on
  its own `(dataset, plate)`) → author dose bins (0–9,
  `assign_pert_dose_corrected`) → the DToxS arms' 48 h exposure encoded at
  v1's 24 h slot (`pert_time 24.0`, `time_idx 2`; gse244740/gse262419 are
  already 24 h). The single-cell encoding is `cell_iname="cardiomyocyte"`,
  `cell_idx=0`, `tissue_idx=0` (num_cell_id 1).
- **Dose sources per dataset**: dtoxs_primary `primary_concentrations.py`
  (paper Table-1 KIs, GSE217423 series-matrix shared non-KIs, literature Cmax
  primary-only drugs; ADR 0015); dtoxs_ipsc its series-matrix concentrations;
  gse244740 the 8-point half-log series from the paper Table-1 tops (already
  in `_datasets.py`); gse262419 real µM from `GSE262419_hash.csv`. **No
  dataset is single-dose any more**: every treated sample carries a positive
  real concentration.
- The union h5ad is schema-identical to v1's `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` in every field
  the XPert loader consumes (obs `pert_id, pert_idx, pert_dose[_idx],
  pert_time[_idx], cell_iname, cell_idx, tissue_idx` + additive `dataset_id`;
  var = 978 landmarks in L1000 `gene_idx` order with `gene_id`/`gene_symbol`;
  `obsm["X_ctl"]` per-sample matched control; X float32 log2(CPM+1) within the
  pretrain bin range `[0, 15.002]`).
- **Delta report** (`output/finetune_delta.md`) audits v3 against v2's union
  (well-identity matched for the DToxS arms): the composition change and the
  per-plate normalization carried over. v1's frozen `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` stays
  the dataset-wide historical anchor.
- Only samples whose compound resolves to an HG drug (`drug_idx ∈ [0, 8981)`)
  enter; the manifest classifies every compound (in_hg / extension / unmatched)
  so exclusions are auditable.

## Finetune-data composition report (step 03)

`03_cm_finetune_data/finetune_report.py` emits a self-contained HTML report
(`output/finetune_report.html`): a **collapsible drug list** grouped by the
author L1000 `moa` (from `pert_model/xpert/metadata/l1000_sdst_drugs_info_8432.csv`,
`Unclassified` for compounds without a recorded moa) and 2D **PCA + UMAP of
cell embeddings** — L1000 pretrain aggregated per cell line (mean embedding) +
every union sample, colored by cell type with the CM samples split by
`dataset_id`. The embedding is the pretrained XPert `cell_Embeddings` output:
binned expression (128 bins over [0, 15.002]) + per-gene PPI vectors,
LayerNorm, mean-pooled over the 978 landmarks.

## Key invariants (gates)

- **Step 01 parity** — `01_drug_embed/output/all_drugs_idx2HG.npy` is
  byte-identical to the author reference (`np.array_equal`); `drug_list.csv` is
  1,004 rows (identical to v1).
- **Step 02 parity** — `02_cm_baseline/output/cm_ctl_landmark.csv` reproduces
  v1's DToxS-primary control consensus.
- **Union schema gate** — `03_cm_finetune_data/schema_gate.py` asserts the
  union h5ad is structure-identical to v1's `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` in the fields
  the XPert loader consumes, **and** the data-semantics invariants (author time
  `{3,6,24}` / `time_idx 2`, **no zero-dose treated samples** (real
  concentrations), X within the pretrain range `[0, 15.002]`, union composition
  = the four finetune datasets, manifest `n_union` matches union rows).
- **Per-plate gates (ADR 0016)** — every treated plate carries ≥1 control well
  on it (0 samples excluded); each treated sample's `X_ctl` equals the mean of
  its own `(dataset, plate)` control rows (asserted in step 03).
- **Fold isolation gates** — `05_xpert_finetune_bench/create_cm_splits.py`:
  1) primary test sets **byte-identical to v1.1's** for all 20 fold columns
  (matched by well id), 2) no identical well across train+test per fold,
  3) cold-drug test drugs absent from **all four datasets'** train,
  4) gse244740/gse262419 rows are train-only (never `test`). Plus the ADR 0014
  invariant checks over the **benchmark surface** (primary + iPSC rows):
  value set subset of {train,test,''}, warm-start overlap 1.0, per-fold test
  fraction in [0.10, 0.30], every benchmark drug and every primary sample
  covered, no stale columns.
- **Gene gates** — `_lib/verify_genes.py` per dataset: G1 all present landmarks
  resolve (gse262419 platform-missing = {WRB, FAM57A} only), G2 all legacy
  aliases author-metadata verified, G3 union-column audit (finetune datasets
  only). Any G1/G2 failure aborts (exit 1).
- **Drug-inclusion gate** — every union obs `pert_idx` lies in `[0, 8981)`;
  the manifest is complete and exclusions match the h5ad row counts.
- Steps 05/06/07 train on the union h5ad only; the frozen v1 `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`
  is used solely as the schema reference; v1.1's `dtoxs_lincs.h5ad` is the
  fold reference; v2's `cm_union.h5ad` is the delta reference (all read-only).

## Entry points

```bash
cd workflow/xpert-fda-miri-v3
python 01_drug_embed/drug_list.py                # 1,004 L1000-contained FDA drugs
python 01_drug_embed/drug_embeddings.py          # author-reference HG (no HG training)
python _lib/verify_genes.py                      # per-dataset gene gates (G1/G2 anytime; G3 needs step-03 union h5ad)
python 02_cm_baseline/cm_baseline.py             # DToxS-primary consensus baseline
python 03_cm_finetune_data/build_union.py        # per-dataset h5ads + union + manifest + delta + schema gate
python 03_cm_finetune_data/finetune_report.py    # composition report (HTML)
python 05_xpert_finetune_bench/create_cm_splits.py  # import v1.1 folds + iPSC/gse244740/gse262419 isolation + gates
# GPU steps (05 benchmark / 06 finetune / 07 predict) then 08–11
python 11_integrate/integrate.py
python _lib/build_screen_notebook.py             # regenerate screen.ipynb
```

## Gotchas

- Never regenerate UniMol features from SMILES at predict time — step 07
  indexes the precomputed array by `drug_idx`.
- The HG table is 8981×256 (author reference only); do not point the config at
  an 11104-row extension table.
- The finetune/benchmark config (`06_xpert_finetune/config/cm.yaml`) uses the
  `cm` dataset block and CLI dataset name `cm`; the XPert library treats `cm`
  like `dtoxs` (additive branches in `utils.py`, `model_xpert.py`,
  `predict.py`). Existing `dtoxs` behavior is unchanged. The pretrain eval
  config (`05_xpert_finetune_bench/config/cm_pretrain_eval.yaml`) mirrors the
  pretrain dims under the `cm` block (num_pert_dose 10, num_cell_id 240,
  num_tissue_id 21) so the pretrain checkpoint loads tolerantly.
- `create_cm_splits.py` **imports** v1.1's fold columns for primary rows and
  never regenerates them — re-running the splitter on the union would drift the
  test sets. It requires the v1.1 DToxS-LINCS h5ad (`V11_DTOXS_LINCS_H5AD`,
  override with `--v11`) to carry the fold columns.
- Step 07 predicts at `--dose-idx 0 --time-idx 2` (the v1 exposure; the
  library `predict.py` defaults to time 2). After any step-05/06 retrain,
  steps 08–11 must be re-run on the regenerated step-07 DEG table.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`) and `.env` are gitignored
  — never commit them.
- Steps 05/06 re-run the XPert finetune on the union — GPU-fp-nondeterministic
  across environments; v3's archived predictions are environment-specific,
  exactly as in v1/v2.
- v3 has **no report stage** (no `12_report/`); the workflow ends at
  `11_integrate`. v1 archives its report as a non-numbered `report/` stage.
