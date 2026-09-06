# AGENTS.md

The **union-CM finetune MIRI screening workflow** (`workflow/xpert-fda-miri-v2/`;
see docs/adr/0009 + docs/adr/0015): steps 01–11 form the 11-step sibling of the finalized
`workflow/xpert-fda-miri-v1/` (ADR 0007) whose **only substantive change is the
XPert finetune data**. Instead of v1's frozen DToxS-only `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`,
steps 05/06/07 train on the workflow-built **union h5ad** of the primary + iPSC
CM datasets (`dtoxs_primary` from the LINCS LDS-1593 well-level package +
`dtoxs_ipsc` from GSE217423, step 03), **per-plate normalized** since ADR 0016
(matched control + QN per `(dataset, plate)`). Since ADR 0015 the finetune
carries **real-concentration dose bins** (not v1's single-dose convention) with
a single-cell encoding, and the benchmark test sets are **byte-identical to
v1.1's** (fold columns imported verbatim for the primary arm). Everything else —
the screened drug scope, step-02 baseline, MIRI phenotype/scoring, filters, and
integration — is v1-identical, preserved by renumbering the v1 steps to make
room for the prep step. Step 12 (`12_drug_survey/`) is an **agentic evidence
survey of the step-11 top-20 candidates** (one biomni agent per drug) that
writes per-drug bilingual dossiers (`drugs/{drug}.md`) and the aggregate
`reports@20.html` — it is read-only on steps 01–11 and does not change any
screen output.

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
  union h5ad + maps + manifest are regenerated workflow artifacts. (The DToxS
  primary finetune transcriptomes come from the released LINCS LDS-1593
  well-level package `pert_data/dtoxs_primary/Data/Level3`, not the raw GEO
  deposit which was deleted 2026-08-18.)
- The **"no external database" provenance rule is amended for the dose field
  only** (ADR 0015): the primary arm's finetune concentrations are cited
  (paper Table 1 KIs / GSE217423 series-matrix shared non-KIs / literature
  clinical Cmax for the 8 primary-only drugs) in
  `03_cm_finetune_data/primary_concentrations.py`.

## Structure (12-step layout)

- `01_drug_embed/` — merged drug list (`drug_list.py`, 1,004 drugs) + author-reference HG embeddings (`drug_embeddings.py`, byte-parity gate). Identical to v1.
- `02_cm_baseline/` — consensus baseline from the DToxS-primary control wells of the released LINCS LDS-1593 package (`pert_data/dtoxs_primary/Data/Level3`); byte-identical to v1's baseline.
- `03_cm_finetune_data/` — **the v2 change**: builds one finetune h5ad per dataset (`output/dtoxs_primary.h5ad` 662 rows, `output/dtoxs_ipsc.h5ad` 296 rows) and concatenates them into the union `output/cm_union.h5ad` (958 rows), from the two finetune datasets (read-only), **per-plate normalized** since ADR 0016 (QN + matched control per `(dataset, plate)`), plus per-dataset gene maps (`output/gene_maps/`), the compound→drug-index manifest (`output/drug_manifest.csv`), the primary-concentration map (`primary_concentrations.py`, cited), a schema gate vs v1's `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`, the delta report vs v1.1's DToxS-LINCS h5ad (`output/finetune_delta.md`, incl. the per-plate deviation), and the composition report `output/finetune_report.html`.
- `04_xpert_pretrain/` — the workflow-owned pretrain unit (weight + `config/l1000.yaml`), carried over untouched.
- `05_xpert_finetune_bench/` — **two independent benchmark methods** (ADR 0018/0019):
  - **full/** (the official fold benchmark, the ADR 0018 *dtoxs-primary bench* — relocated verbatim under `full/` by ADR 0019, plot kept `full/benchmark.png`): `full/create_cm_splits.py` imports v1.1's `warm-start_fold-1..10` + `cold-drug_fold-1..10` verbatim for primary rows (matched by well id, byte-identical test sets) and assigns iPSC rows under the isolation rules (warm-start all-train; cold-drug test drugs excluded from both datasets' train) + `cm_full` (all-train), with three isolation gates; `full/xpert_finetune_bench.py` benches the 20 folds at **4 epochs** and — as the **first column** (`pretrain`, amber) — zero-shot-evaluates the step-04 L1000 pretrain weight on all 20 fold test sets (no finetuning) under the workflow-owned `config/cm_pretrain_eval.yaml` (step-top-level, shared). Comparable column-by-column with v1.1 (identical test sets). Cached fold weights/results under `full/output/`.
  - **data-ablation/** (ADR 0019, decoupled from the official folds): `xpert_data_ablation.py` measures how the **finetune-corpus scale** affects XPert test performance on v2's own union (958 rows / 45 drugs). Because a row-uniform test would draw ~40 of 45 drugs into test, the per-seed test sets are **drug-partitioned**: a disjoint tripartition of the 45 union drugs into 3×15 test drugs — each seed tests one third, so every union drug is tested exactly once — with a ~15%-of-union unified test set drawn only from that third. 3 seeds × 2 modes (warm-start | cold-drug) × 5 train-size fractions **{20%, 40%, 60%, 80%, 100%}** (warm pool ~814 rows all 45 drugs; cold pool the 30 hold-out drugs, ~610–655 rows) = 30 runs + 3 L1000-pretrain zero-shot anchors; recipe identical to the official bench. Outputs under `data-ablation/` (masks, `ablation.csv`, 4×4 `ablation.png`, `ablation_report.md`, per-run weights/profiles/results); drives `src.train` via temp h5ads — `cm_union.h5ad`'s 21 official fold columns stay intact.
- `06_xpert_finetune/` — all-train union finetune (`xpert_finetune.py`, `config/cm.yaml`) -> the final prediction model `output/cm_cm_full/`.
- `07_xpert_predict/` — perturbed-transcriptome prediction + z-score (`perturb_transcriptome.py`, `zscore.py`); predicts at dose_idx 0 / time_idx 2 (the v1 exposure, unchanged by the real-dose finetune).
- `08_miri_phenotype/` — integrated MIRI reference gene sets + rank reference.
- `09_miri_efficacy_scores/` — the three efficacy scorers + z-components + plots.
- `10_filter/` — WCS<0 + >=6/8 critical-pathway filters.
- `11_integrate/` — integrated score over all drugs, ranked candidates, deterministic report.
- `12_drug_survey/` — **agentic top-20 evidence survey** (read-only on steps 01–11): `prompt.md` (the per-drug bilingual agent brief + two-score rubric), `assemble_reports.py` (deterministic dossier → `reports@20.html` builder). Generated outputs: `drugs/{drug}.md` (one bilingual dossier per top-20 drug, each beginning with a parseable key-value header block) + `reports@20.html` (bilingual summary table + embedded dossiers). gitignored (regenerable external-agent output).
- `input/` — frozen assets (author HG reference, UniMol array, idx→SMILES, FDA drug list, curated `MIRI_pathways.json`), carried over from v1.
- `_lib/utils.py` — workflow-root-aware paths; `CM_DATA_ROOT`, `CM_DATASETS` (all five, for the data-lake audit), `FINETUNE_DATASETS = ["dtoxs_primary", "dtoxs_ipsc"]`, `N_L1000_DRUGS = 8981`, `AUTHOR_TIME_MAP = {3:0, 6:1, 24:2}`, `L1000_PRETRAIN_H5AD`, `AUTHOR_HG_REFERENCE`, `V11_DTOXS_LINCS_H5AD` (v1.1 fold/delta source).
- `_lib/gene_mapping.py` — **per-dataset** CM→L1000 landmark maps (symbols / Ensembl / SYMBOL_entrezid), legacy-symbol aware, `corrected_gene_name` never a match key.
- `_lib/drug_resolution.py` — compound→L1000 drug index (canonical-SMILES match, else author `cmap_name` name match) + the auditable manifest.
- `_lib/verify_genes.py` — the extended per-dataset gene gates (G1 completeness, G2 alias verification; G3 union-column audit for the finetune datasets only), writes `_lib/gene_mapping_report.md`.

## The union finetune data (step 03, ADR 0015)

- The finetune scope is **primary + iPSC only** (`dtoxs_primary` from the LINCS
  LDS-1593 well-level package, `dtoxs_ipsc` from GSE217423); gse244740 /
  gse262419 / gse247541 are **no longer extracted in the build** — the
  standalone drug-overlap report (`03_cm_finetune_data/drug_overlap_report.py`)
  still audits all five for the data lake. Its summary table compares, per
  dataset, treated 化合物总数 / 有 SMILES 数 / in-HG 药物数 (+ treated 样本数 and a
  5-dataset union row), and its gse262419 well-join is **plate-aware**
  (`(Plate_ID, well)` keyed against `GSE262419_hash.csv`) since 2026-08-27 —
  the old well-only join collapsed the 464 chemicals to 33 names and
  undercounted gse262419's in-HG drugs (26 → 101).
- **Per-plate normalization (ADR 0016)** — the QN alignment and the matched
  control are grouped per `(dataset_id, plate)` (one monotonic transform per
  plate, before matched-control computation): each sample's control baseline is
  the mean of the control wells on its own plate (dtoxs_primary CTRL /
  CTRL-DMSO / CTRL-EtOH on plates 1/2/4/5; dtoxs_ipsc CTRL on plate 6). Every
  treated plate has ≥1 control well (0 samples excluded), asserted in step 03.
  v1's frozen `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` stays dataset-wide (the historical anchor);
  the deviation is recorded in the delta report.
- Each dataset is preprocessed independently: 978 landmarks → CPM → log2(CPM+1)
  → per-`(dataset_id, plate)` quantile alignment to the L1000 pretrain pooled
  distribution (QN reference cached at `output/qn_reference.npy`) → per-plate
  matched control → author dose bins (0–9, `assign_pert_dose_corrected`) → the
  48 h exposure encoded at v1's 24 h slot (`pert_time 24.0`, `time_idx 2`).
  The single-cell encoding is `cell_iname="cardiomyocyte"`, `cell_idx=0`,
  `tissue_idx=0` (num_cell_id 1; the model never distinguishes donors/lines).
- **Real concentrations** (ADR 0015): the primary arm's dose comes from
  `primary_concentrations.py` — paper Table 1 clinical concentrations for the
  23 KIs (SOR 0.5 uM), the GSE217423 series matrix for the 15 shared non-KIs,
  and literature-curated clinical Cmax with citations for the 8 primary-only
  drugs (SB-431542 = in-vitro 10 uM). The iPSC arm keeps its series-matrix
  concentrations (SOR 1.0 uM — the primary 0.5 vs iPSC 1.0 SOR divergence is
  intended and documented). **No dataset is single-dose any more**: every
  treated sample carries a positive real concentration, so the old
  "zero dose only in dtoxs_primary" schema assertion is gone.
- The union h5ad is schema-identical to v1's `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` in every field
  the XPert loader consumes (obs `pert_id, pert_idx, pert_dose[_idx],
  pert_time[_idx], cell_iname, cell_idx, tissue_idx` + additive `dataset_id`;
  var = 978 landmarks in L1000 `gene_idx` order with `gene_id`/`gene_symbol`;
  `obsm["X_ctl"]` per-sample matched control; X float32 log2(CPM+1) within the
  pretrain bin range `[0, 15.002]`). Union = 662 + 296 = **958 rows**,
  45 unique drug_idx (46 primary compounds; doxorubicin/epirubicin share 4556).
- **Delta report** (`output/finetune_delta.md`) audits the primary arm against
  v1.1's `dtoxs_lincs.h5ad` (well-identity matched, 662/662): the dose_idx
  re-binning (all bin 0 → the real-concentration distribution) and the
  `cell_iname` change (`cardio` → `cardiomyocyte`). The iPSC arm is new in v2.
- Only samples whose compound resolves to an HG drug (`drug_idx ∈ [0, 8981)`)
  enter; the manifest classifies every compound (in_hg / extension / unmatched)
  so exclusions are auditable (primary non-mapped treated codes such as TRS /
  DIC / END / SUN+ / TNF / IGF are counted as `n_excl_nodose`).

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
  = the two finetune datasets, manifest `n_union` matches union rows).
- **Per-plate gates (ADR 0016)** — step 03: every treated plate carries ≥1
  control well on it (0 samples excluded); each treated sample's `X_ctl`
  equals the mean of its own `(dataset, plate)` control rows.
- **Fold isolation gates** — `05_xpert_finetune_bench/full/create_cm_splits.py`:
  1) primary test sets **byte-identical to v1.1's** for all 20 fold columns
  (matched by well id), 2) no identical well across train+test per fold,
  3) cold-drug test drugs absent from both datasets' train. Plus the ADR 0014
  invariant checks (value set, warm-start overlap 1.0, per-fold test fraction
  in [0.10, 0.30], every primary sample/drug covered, no stale columns).
- **Gene gates** — `_lib/verify_genes.py` per dataset: G1 all present landmarks
  resolve (gse262419 platform-missing = {WRB, FAM57A} only), G2 all legacy
  aliases author-metadata verified, G3 union-column audit (finetune datasets
  only). Any G1/G2 failure aborts (exit 1).
- **Drug-inclusion gate** — every union obs `pert_idx` lies in `[0, 8981)`;
  the manifest is complete and exclusions match the h5ad row counts.
- Steps 05/06/07 train on the union h5ad only; the frozen v1 `workflow/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`
  is used solely as the schema reference; v1.1's `dtoxs_lincs.h5ad` is the
  fold/delta reference (read-only).

## Entry points

```bash
cd workflow/xpert-fda-miri-v2
python 01_drug_embed/drug_list.py                # 1,004 L1000-contained FDA drugs
python 01_drug_embed/drug_embeddings.py          # author-reference HG (no HG training)
python _lib/verify_genes.py                      # per-dataset gene gates (G1/G2 anytime; G3 needs step-03 union h5ad)
python 02_cm_baseline/cm_baseline.py             # DToxS-primary consensus baseline
python 03_cm_finetune_data/build_union.py        # per-dataset h5ads + union + manifest + delta + schema gate (per-plate, ADR 0016)
python 03_cm_finetune_data/finetune_report.py    # composition report (HTML)
python 05_xpert_finetune_bench/full/create_cm_splits.py  # import v1.1 folds + iPSC isolation + gates
python 05_xpert_finetune_bench/full/xpert_finetune_bench.py  # official fold bench -> full/ (GPU)
python 05_xpert_finetune_bench/data-ablation/xpert_data_ablation.py --dry-run  # drug-partitioned masks + isolation gates (no GPU)
python 05_xpert_finetune_bench/data-ablation/xpert_data_ablation.py  # data-ablation bench -> data-ablation/ (GPU, 30 runs)
# GPU steps (05 benchmark / 06 finetune / 07 predict) then 08–11
python 11_integrate/integrate.py
# Step 12 (agentic, external biomni; read-only on steps 01–11):
#   1) run one @biomni subagent per top-20 drug per prompt.md  -> drugs/{drug}.md
#   2) python 12_drug_survey/assemble_reports.py              -> reports@20.html
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
  config (`05_xpert_finetune_bench/config/cm_pretrain_eval.yaml`) mirrors
  the pretrain dims under the `cm` block (num_pert_dose 10, num_cell_id 240,
  num_tissue_id 21) so the pretrain checkpoint loads tolerantly. The
  data-ablation bench shares it (no copy) and generates its own per-run YAMLs
  under `05_xpert_finetune_bench/config/` (the same nesting as the fold bench's).
- `full/create_cm_splits.py` **imports** v1.1's fold columns for primary rows
  and never regenerates them — re-running the splitter on the union would
  drift the test sets. It requires the v1.1 DToxS-LINCS h5ad
  (`V11_DTOXS_LINCS_H5AD`, override with `--v11`) to carry the fold columns.
- **Two benchmarks, different test sets** (ADR 0019): the official `full/`
  bench uses the ADR 0014/0015 folds (test sets primary-only, byte-identical
  to v1.1); the `data-ablation/` bench uses its own per-seed **drug-partitioned
  unified test sets** (a disjoint tripartition of the 45 union drugs — every
  drug tested exactly once across the 3 seeds) and is deliberately decoupled.
  Never mix the two — the ablation's 100% point re-runs on its own split
  scheme and is NOT the official full/ bench. The ablation drives `src.train`
  via temporary h5ads and never writes to `cm_union.h5ad` (the 21 official
  fold columns stay intact); reruns are idempotent per (seed, mode, fraction)
  via `data-ablation/` caches, resumable after interruption.
- `xpert_data_ablation.py` is a **diagnostic analysis** only: it does not
  change the final screen model (step 06 `cm_full`) or any downstream step.
- Step 07 predicts at `--dose-idx 0 --time-idx 2` (the v1 exposure; the
  library `predict.py` defaults to time 2). After any step-05/06 retrain,
  steps 08–11 must be re-run on the regenerated step-07 DEG table.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`) and `.env` are gitignored
  — never commit them.
- Steps 05/06 re-run the XPert finetune on the union — GPU-fp-nondeterministic
  across environments; v2's archived predictions are environment-specific,
  exactly as in v1.
- Step 12 (`12_drug_survey/`) is **agentic and external**: one biomni agent
  per top-20 drug (prompt.md rubric; agent proposes `score_miri`/`score_iri`,
  the step owner confirms against the cited evidence). Never fabricate a report
  or citation — agents report "no evidence found" honestly. Generated dossiers
  + `reports@20.html` are gitignored (regenerable external-agent output) and
  resumable (an existing `drugs/{drug}.md` is not re-researched). Step 12 does
  not alter any step-01–11 output; `assemble_reports.py` flags missing/invalid
  dossier headers rather than silently dropping rows.
