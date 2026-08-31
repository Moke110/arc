# AGENTS.md

The **finalized MIRI screening workflow** (`workflow/xpert-fda-miri-v1/`, formerly `xpert-fda-v3`; see docs/adr/0007): an **author-aligned, L1000-contained FDA screen**. It screens **only the FDA drugs already present in the L1000 embedding table** (exact canonical-SMILES match, `drug_idx < 8981`; 1,004 of 3,127 drugs). Because every screened drug uses the author-reference HG embedding directly, **no HG training / low-LR finetune / FDA-extension inference is ever run** — the workflow is fully aligned with the author's original L1000-based XPert.

## Scope and provenance (author + DToxS only)

- The v1/v2 FDA-extension variants (2,123 FDA drugs via HG finetune / author-pinned inference) are **removed**; the orphaned FDA-extension data/weights live in the gitignored archive `reference/hg-fda-extension-archive/`.
- **Drug info** comes only from the workflow-owned `input/` files (author L1000 SMILES/HG tables, precomputed UniMol array, inherited FDA query list `fda_drugs.csv`). The workflow never queries an external drug database and never generates drug features from SMILES.
- **Gene info** comes only from the XPert-author `pert_model/xpert/metadata/l1000_gene_info_978.csv` and the DToxS package. No external gene databases.
- The 1,004 screened drugs are a strict subset of the author's L1000 table; their embeddings are the author reference byte-identical (step 01 parity gate).
- The frozen DToxS finetune h5ad (`03_cm_finetune_data/output/dtoxs.h5ad`) is workflow-owned (moved from `pert_data/dtoxs_primary/v1-frozen/`); the well-level package stays frozen in `pert_data` for the step-02 baseline.

## Structure (finalized 11-step layout, identical to v1.1/v2/v3)

- `01_drug_embed/` — merged drug list (`drug_list.py`) + author-reference HG embeddings (`drug_embeddings.py`, byte-parity gate).
- `02_cm_baseline/` — consensus cardiomyocyte baseline (`cm_baseline.py`), authoritative gene mapping.
- `03_cm_finetune_data/` — owns the frozen DToxS finetune data `output/dtoxs.h5ad` (521 samples / 36 drugs), moved here from `pert_data/dtoxs_primary/v1-frozen/` so the workflow is self-contained on its finetune data (matching the v1.1/v2/v3 step-03 data prep). Not rebuilt; its fold columns are regenerated in place by step 05.
- `04_xpert_pretrain/` — the workflow-owned pretrain unit: used pretrain weight + `config/l1000.yaml` (data roots point at `pert_data/l1000-xpert/data/` and `input/`).
- `05_xpert_finetune_bench/` — XPert DTOXS-split benchmark on the **unified folds** (ADR 0014): `create_dtoxs_splits.py` regenerates `warm-start_fold-1..10` (drug-overlap) + `cold-drug_fold-1..10` (drug-disjoint) + `dtoxs_full` on the frozen `03_cm_finetune_data/output/dtoxs.h5ad` in place (invariant gate, seed 2024); `xpert_finetune_bench.py` benches the 20 folds and overlays per-fold dots on the mean ± std bars, and — as the **first column** (`pretrain`, amber) — zero-shot-evaluates the step-04 L1000 pretrain weight on all 20 fold test sets (no DToxS finetune) under the workflow-owned `config/dtoxs_pretrain_eval.yaml` (pretrain architecture/bin `[0, 15.002]`; the drifted aux `class_fc` head is dropped by the library's tolerant test-mode load).
- `06_xpert_finetune/` — all-train DTOXS finetune (`xpert_finetune.py`) + step-owned `config/dtoxs.yaml`.
- `07_xpert_predict/` — merged perturbed-transcriptome prediction (`perturb_transcriptome.py`) + z-score (`zscore.py`).
- `08_miri_phenotype/` — integrated MIRI reference gene sets + rank reference (`miri_phenotype.py`).
- `09_miri_efficacy_scores/` — the three efficacy scorers + z-components + plots (`miri_efficacy_scores.py`).
- `10_filter/` — WCS<0 + >=6/8 critical-pathway filters (`filter.py`).
- `11_integrate/` — integrated score over all drugs, ranked candidates, deterministic report (`integrate.py`).
- `report/` — **archived, non-numbered stage**: agentic bilingual final drug report (`run_drug_agents.py`, `assemble_report.py`, `assemble_cn_report.py`, `prompt.md`; per-drug JSON cache in `cache/`, outputs in `output/`). Kept out of the numbered steps so v1's layout is identical to v1.1/v2/v3.
- `input/` — **real directory** owning the frozen UniMol array (`all_drugs_unimol_arr.npy`), idx→SMILES mapping (`all_drugs_idx2smi.npy`), FDA drug list, curated `MIRI_pathways.json`, and the **author HG reference** (`HG_drug_embeddings.npy`). The DToxS finetune h5ad now lives in step 03.
- `_lib/utils.py` — workflow-root-aware paths; `N_L1000_DRUGS = 8981`; `resolve_l1000_drug_rows()`; `AUTHOR_HG_REFERENCE = input/HG_drug_embeddings.npy`.
- `_lib/scoring.py` — the three efficacy scorers + integration.
- `_lib/gene_mapping.py` — **authoritative DToxS→L1000 landmark gene map**, derived from **XPert-author + DToxS data only**. Resolves the 18 legacy-symbol landmarks (e.g. WASHC5 ← DToxS KIAA0196, ELP1 ← IKBKAP, JPT2 ← HN1L) using the author's own `l1000_gene_info_978.csv` `gene_name` column. `corrected_gene_name` is never a match key — it rejects the JPT2 ← CRAMP1 trap (CRAMP1 is Entrez 57585, a *different* gene; the DToxS table carries HN1L, JPT2's real alias, Entrez 90861).
- `_lib/verify_genes.py` — the gene-mapping gate (G1 completeness, G2 author-metadata alias verification, G3 frozen-h5ad audit), writes `_lib/gene_mapping_report.md`.

## Gene matching (why it matters)

DToxS (2015-era 3'-DGE symbols) contains all 978 L1000 landmark genes, but 18 landmarks only match via a **retired/alternative symbol**. Direct symbol matching zeroes them (v1's step-02 baseline lost 29 genes). The finalized baseline maps **978/978** landmarks and only the 11 genuinely undetected genes (SKIV2L, MIF, UBR7, PAK6, HSPA1A, BTK, PSMB8, PLEKHM1, HLA-DMA, PLSCR3, TSTA3) stay at 0.

## Key invariants (gates)

- **Step 01 parity** — `01_drug_embed/output/all_drugs_idx2HG.npy` is byte-identical to the author reference (`np.array_equal`).
- **Step 01 scope** — every drug in `01_drug_embed/output/drug_list.csv` has `drug_idx < 8981` (exact SMILES match). No FDA-extension rows are ever predicted.
- **Gene gates** — `_lib/verify_genes.py`: G1 all 978 landmarks resolve; G2 all 18 aliases author-metadata verified; G3 no confirmed h5ad gene misassignment. Any G1/G2 failure aborts (exit 1); G3 findings are reported.
- **Fold gate** — `create_dtoxs_splits.py` asserts after writing: warm-start drug overlap 1.0, cold-drug drug overlap 0.0, per-fold test fraction ∈ [0.10, 0.30], full sample/drug coverage, no stale fold columns (ADR 0014).
- The DToxS finetune (steps 05/06) uses the workflow-owned frozen `03_cm_finetune_data/output/dtoxs.h5ad` (not rebuilt; its fold columns are regenerated in place under ADR 0014, re-basing the archived benchmark).

## Entry points

```bash
cd workflow/xpert-fda-miri-v1
python 01_drug_embed/drug_list.py                # 1,004 L1000-contained FDA drugs
python 01_drug_embed/drug_embeddings.py          # author-reference HG (no HG training)
python _lib/verify_genes.py                      # gene-mapping gates (run before step 02)
python 02_cm_baseline/cm_baseline.py
python 05_xpert_finetune_bench/create_dtoxs_splits.py  # unified folds + invariant gate
# GPU steps (05 benchmark / 06 finetune / 07 predict) then 08–11
python 11_integrate/integrate.py
# archived agentic report (external biomni backend) — resumable
python report/run_drug_agents.py
python report/assemble_report.py
python report/assemble_cn_report.py
python _lib/build_screen_notebook.py             # regenerate screen.ipynb
```

## Gotchas

- Never regenerate UniMol features from SMILES at predict time — step 07 indexes the precomputed array by `drug_idx`.
- The HG table is 8981×256 (author reference only); do not point the config at an 11104-row extension table.
- Pretrain weight + both YAMLs (pretrain + finetune) are workflow-owned; the model library resolves configs against the workflow root.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`), `report/cache/` JSONs, and `.env` are gitignored — never commit them.
- Steps 05/06 re-run the XPert DToxS finetune/benchmark on the author-reference HG — GPU-fp-nondeterministic across environments. The `pretrain` benchmark column re-evaluates the step-04 L1000 weight (no finetune) and is deterministic per GPU; the DToxS-specific embedding tables (dose/time/position) it drops are untrained in that column, exactly as they are at finetune epoch 0.
- The `report/` biomni agent backend (`/home/chang/cc-research`, overridable via `BIOMNI_CC_DIR`) is an external dependency; it is not vendored.
