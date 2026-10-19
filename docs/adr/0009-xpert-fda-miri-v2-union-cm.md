# Union-CM finetune: xpert-fda-miri-v2

The finalized MIRI screen (`stage1/xpert-fda-miri-v1`, ADR 0007) is cloned
into a **union-CM finetune sibling** `stage1/xpert-fda-miri-v2` whose only
substantive change is the XPert finetune data: instead of the frozen
DToxS-only `input/dtoxs.h5ad`, the finetune/benchmark/predict steps train on a
workflow-built **union h5ad** of all five `pert_data/cm` datasets, produced by
a new step (`03_cm_finetune_data`). The rest of the screen — drug scope,
baseline, MIRI scoring, filters, integration, report — is v1-identical.

## Context

v1 finetunes XPert on a single DToxS finetune dataset (521 samples, 26 FDA
kinase inhibitors) while five cardiomyocyte (CM) drug-perturbation datasets
(`dtoxs_primary`, `dtoxs_ipsc`, `gse244740`, `gse247541`, `gse262419`) now
live together under `pert_data/cm/` — spanning primary CMs, iPSC-CMs and
hiPSC-CMs with ~40× more perturbed samples, far more compounds, and multiple
dose/time points. The researcher wants the drug-response model that powers the
FDA screen trained on the full CM transcriptomic corpus, not just DToxS.

## Decision

- **Workflow scaffold** — `stage1/xpert-fda-miri-v2` is a copy of v1 with a
  new prep step inserted, giving a 12-step layout:
  `01_drug_embed`, `02_cm_baseline`, **`03_cm_finetune_data`** (new),
  `04_xpert_pretrain`, `05_xpert_finetune_bench`, `06_xpert_finetune`,
  `07_xpert_predict`, `08_miri_phenotype`, `09_miri_efficacy_scores`,
  `10_filter`, `11_integrate`, `12_report`. Only numeric prefixes shift; the
  descriptive steps are carried over unchanged (v1 parity by construction).
- **Union finetune data (step 03)** — reads `pert_data/cm` **read-only** and
  emits `03_cm_finetune_data/output/cm_union.h5ad`, schema-identical to v1's
  `input/dtoxs.h5ad` in every field the XPert loader consumes (verified by the
  step's schema gate). Only samples whose compound resolves to an HG drug
  (`drug_idx ∈ [0, 8981)`) enter; a per-compound manifest classifies every
  dataset compound (in_hg / extension / unmatched) so exclusions are auditable.
- **Drug resolution** — exact canonical-SMILES match (connectivity only,
  `isomericSmiles=False`, matching v1's indices) where the dataset or DToxS
  metadata provides a SMILES, else drug-name match against the author L1000
  `cmap_name` table (normalized; salt/racemic-prefix variants handled as pure
  name normalization). No external database, no SMILES regeneration.
- **Gene resolution** — per-dataset landmark maps in each dataset's own
  identifier space (symbols / Ensembl / SYMBOL_entrezid) with legacy-symbol
  awareness from the author `l1000_gene_info_978.csv`; `corrected_gene_name`
  is never a match key (the JPT2←CRAMP1 trap). `gse262419`'s TempO-Seq panel
  provably lacks WRB/FAM57A (probe rows point at different Entrez genes);
  they are documented platform-missing and zeroed, not silently mismatched.
- **Dose/time/cell/tissue indices** — precomputed in step 03: dose via the
  XPert `assign_pert_dose_corrected` bins (0–9), time via a step-owned union
  map `{24:0, 48:1, 144:2}` (the library fallback only covers {3,6,24} h), all
  CM → `cell_idx 0 / tissue_idx 0` (`num_cell_id 1 / num_tissue_id 1`).
  The screen predicts at `time_idx 0` (= 24 h, the L1000-standard exposure v1
  conditioned on) so only the finetune data differs between v1 and v2.
- **Config** — the finetune/benchmark YAML uses a `cm` dataset block
  (`num_pert_dose 10`, `num_pert_time 3`, `num_cell_id 1`, `num_tissue_id 1`,
  `max_value`/`min_value` from the union distribution) and CLI dataset name
  `cm`. This requires a **minimal additive XPert-library extension** — treat
  `cm` like `dtoxs` in `_load_data_and_drug_feat` (utils.py), the dose/time
  embedding branch (model_xpert.py), and the config-block selection in
  predict.py. Existing `dtoxs` behavior is unchanged; this is a deliberate,
  documented deviation from v2's "no library changes" scope, forced by the
  ticket-05 requirement that no `dtoxs` reference remains in the config.
- **Benchmark scope** — the step-05 benchmark runs the drug-stratified
  **warm-start** folds (seed 2024) on the union; the **`cm_full`** all-train
  split is declared in the config and trained by step 06 (an all-train split
  has no held-out test set, so it is not evaluated by the benchmark).
  `cold-drug`/`cold-cate` are dropped: drug categories exist only for the
  DToxS KI panel, and the union's compounds are heterogeneous. The final
  prediction model (step 06) trains on `cm_full`.
- **Unchanged steps are copies** — steps 01, 02, and 08–12 (renumbered) are
  carried over with only path/step-number updates. Step-02's baseline is
  rebuilt from the DToxS-primary raw GEO deposit (the old `pert_data/dtoxs`
  path no longer exists) and reproduces v1's cached baseline byte-for-byte.
- **Docs** — a v2 `AGENTS.md`, updated root `AGENTS.md`/`CONTEXT.md`, and this
  ADR record the v1/v2 relationship.

## Consequences

- Two named MIRI workflows coexist: v1 (DToxS-only finetune) and v2
  (union-CM finetune); the rest of the screen is identical, so v1/v2 output
  differences isolate the finetune-data effect.
- The union is platform- and batch-heterogeneous; `gse262419` alone is ~70%
  of the samples. Per-dataset `log2(CPM+1)` normalization does not remove
  batch/platform shifts; platform-effect mitigation is a future decision.
- The union h5ad, gene maps, drug manifest, and schema gate are regenerated
  workflow artifacts (gitignored); `pert_data/cm` is never modified.
- GPU steps (05/06/07) are GPU-fp-nondeterministic across environments; v2's
  archived predictions are environment-specific, exactly as in v1.
