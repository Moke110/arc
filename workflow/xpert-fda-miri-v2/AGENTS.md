# AGENTS.md

The **union-CM finetune MIRI screening workflow** (`workflow/xpert-fda-miri-v2/`;
see docs/adr/0009): the 12-step sibling of the finalized
`workflow/xpert-fda-miri-v1/` (ADR 0007) whose **only substantive change is the
XPert finetune data**. Instead of v1's frozen DToxS-only `input/dtoxs.h5ad`,
steps 05/06/07 train on the workflow-built **union h5ad** of all five
`pert_data/cm` datasets (step 03). Everything else — the screened drug scope,
step-02 baseline, MIRI phenotype/scoring, filters, integration, and the
agentic report — is v1-identical, preserved by renumbering the v1 steps to make
room for the prep step.

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
- `pert_data/cm/` is **read-only**: step 03 never modifies a source file; the
  union h5ad + maps + manifest are regenerated workflow artifacts.

## Structure (12-step layout)

- `01_drug_embed/` — merged drug list (`drug_list.py`, 1,004 drugs) + author-reference HG embeddings (`drug_embeddings.py`, byte-parity gate). Identical to v1.
- `02_cm_baseline/` — consensus baseline from the DToxS-primary control wells of the raw GEO deposit (`pert_data/cm/dtoxs_primary`); byte-identical to v1's baseline.
- `03_cm_finetune_data/` — **the v2 change**: builds the union finetune h5ad `output/cm_union.h5ad` from all five `pert_data/cm` datasets (read-only), plus per-dataset gene maps (`output/gene_maps/`), the compound→drug-index manifest (`output/drug_manifest.csv`), and a schema gate vs v1's `input/dtoxs.h5ad`.
- `04_xpert_pretrain/` — the workflow-owned pretrain unit (weight + `config/l1000.yaml`), carried over untouched.
- `05_xpert_finetune_bench/` — union warm-start benchmark (`create_cm_splits.py`, `xpert_finetune_bench.py`); `cold-drug`/`cold-cate` dropped (no meaningful union drug categories).
- `06_xpert_finetune/` — all-train union finetune (`xpert_finetune.py`, `config/cm.yaml`) -> the final prediction model `output/cm_cm_full/`.
- `07_xpert_predict/` — perturbed-transcriptome prediction + z-score (`perturb_transcriptome.py`, `zscore.py`).
- `08_miri_phenotype/` — integrated MIRI reference gene sets + rank reference.
- `09_miri_efficacy_scores/` — the three efficacy scorers + z-components + plots.
- `10_filter/` — WCS<0 + >=6/8 critical-pathway filters.
- `11_integrate/` — integrated score over all drugs, ranked candidates, deterministic report.
- `12_report/` — agentic bilingual final drug report (external biomni backend).
- `input/` — frozen assets (author HG reference, UniMol array, idx→SMILES, FDA drug list, curated `MIRI_pathways.json`), carried over from v1.
- `_lib/utils.py` — workflow-root-aware paths; `CM_DATA_ROOT`, `N_L1000_DRUGS = 8981`, `UNION_TIME_MAP = {24:0, 48:1, 144:2}`, `AUTHOR_HG_REFERENCE`.
- `_lib/gene_mapping.py` — **per-dataset** CM→L1000 landmark maps (symbols / Ensembl / SYMBOL_entrezid), legacy-symbol aware, `corrected_gene_name` never a match key.
- `_lib/drug_resolution.py` — compound→L1000 drug index (canonical-SMILES match, else author `cmap_name` name match) + the auditable manifest.
- `_lib/verify_genes.py` — the extended per-dataset gene gates (G1 completeness, G2 alias verification, G3 drug-level audit), writes `_lib/gene_mapping_report.md`.

## The union finetune data (step 03)

- All five `pert_data/cm` datasets enter: `dtoxs_primary` (raw GEO DGE, 48 h),
  `dtoxs_ipsc` (GSE217423 DGE, 48 h), `gse244740` (ScreenSeq, 24 h, 8-point
  half-log dose series), `gse247541` (bulk RNA-seq, 144 h), `gse262419`
  (TempO-Seq, 24 h, ~55% of union samples). gse244740's rows are
  Ensembl-indexed, so it resolves landmarks by Ensembl ID (the spec's
  "symbols" phrasing for gse244740 was a text error).
- Only samples whose compound resolves to an HG drug (`drug_idx ∈ [0, 8981)`)
  enter; the manifest classifies every compound (in_hg / extension / unmatched)
  so exclusions are auditable.
- The union h5ad is schema-identical to v1's `input/dtoxs.h5ad` in every field
  the XPert loader consumes (obs `pert_id, pert_idx, pert_dose[_idx],
  pert_time[_idx], cell_iname, cell_idx, tissue_idx` + additive `dataset_id`;
  var = 978 landmarks in L1000 `gene_idx` order with `gene_id`/`gene_symbol`;
  `obsm["X_ctl"]` per-sample matched control; X float32 log2(CPM+1)).
- Dose/time indices are precomputed in step 03 (not left to the library's
  fixed {3,6,24} h map): the union time map is `{24:0, 48:1, 144:2}` and the
  screen predicts at time_idx 0 (= 24 h, v1's exposure).
- `gse262419`'s TempO-Seq panel provably lacks WRB and FAM57A (their probe rows
  point at different Entrez genes); those two landmarks are zeroed for that
  dataset and documented, not silently mismatched.

## Key invariants (gates)

- **Step 01 parity** — `01_drug_embed/output/all_drugs_idx2HG.npy` is
  byte-identical to the author reference (`np.array_equal`); `drug_list.csv` is
  1,004 rows (identical to v1).
- **Step 02 parity** — `02_cm_baseline/output/cm_ctl_landmark.csv` reproduces
  v1's DToxS-primary control consensus.
- **Union schema gate** — `03_cm_finetune_data/schema_gate.py` asserts the
  union h5ad is structure-identical to v1's `input/dtoxs.h5ad` in the fields
  the XPert loader consumes.
- **Gene gates** — `_lib/verify_genes.py` per dataset: G1 all present landmarks
  resolve (gse262419 platform-missing = {WRB, FAM57A} only), G2 all legacy
  aliases author-metadata verified, G3 no confirmed union-column
  misassignment. Any G1/G2 failure aborts (exit 1).
- **Drug-inclusion gate** — every union obs `pert_idx` lies in `[0, 8981)`;
  the manifest is complete and exclusions match the h5ad row counts.
- Steps 05/06/07 train on the union h5ad only; the frozen v1 `input/dtoxs.h5ad`
  is used solely as the schema reference.

## Entry points

```bash
cd workflow/xpert-fda-miri-v2
python 01_drug_embed/drug_list.py                # 1,004 L1000-contained FDA drugs
python 01_drug_embed/drug_embeddings.py          # author-reference HG (no HG training)
python _lib/verify_genes.py                      # per-dataset gene gates (G1/G2 anytime; G3 needs step-03 union h5ad)
python 02_cm_baseline/cm_baseline.py             # DToxS-primary consensus baseline
python 03_cm_finetune_data/build_union.py        # union h5ad + manifest + schema gate
# GPU steps (05 benchmark / 06 finetune / 07 predict) then 08–12
python 11_integrate/integrate.py
# agentic report (external biomni backend) — resumable
python 12_report/run_drug_agents.py
python 12_report/assemble_report.py
python 12_report/assemble_cn_report.py
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
  `predict.py`). Existing `dtoxs` behavior is unchanged.
- `create_cm_splits.py` normalizes the union obs index to a RangeIndex before
  writing split columns (anndata stores a string index on disk).
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`), `12_report/cache/` JSONs,
  and `.env` are gitignored — never commit them.
- Steps 05/06 re-run the XPert finetune on the union — GPU-fp-nondeterministic
  across environments; v2's archived predictions are environment-specific,
  exactly as in v1.
- The 12_report biomni agent backend (`/home/chang/cc-research`, overridable
  via `BIOMNI_CC_DIR`) is an external dependency; it is not vendored.
