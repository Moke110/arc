# AGENTS.md

The **finalized MIRI screening workflow** (`workflow/xpert-fda-miri-v1/`, formerly `xpert-fda-v3`; see docs/adr/0007): an **author-aligned, L1000-contained FDA screen**. It screens **only the FDA drugs already present in the L1000 embedding table** (exact canonical-SMILES match, `drug_idx < 8981`; 1,004 of 3,127 drugs). Because every screened drug uses the author-reference HG embedding directly, **no HG training / low-LR finetune / FDA-extension inference is ever run** — the workflow is fully aligned with the author's original L1000-based XPert.

## Scope and provenance (author + DToxS only)

- The v1/v2 FDA-extension variants (2,123 FDA drugs via HG finetune / author-pinned inference) are **removed**; the orphaned FDA-extension data/weights live in the gitignored archive `reference/hg-fda-extension-archive/`.
- **Drug info** comes only from the workflow-owned `input/` files (author L1000 SMILES/HG tables, precomputed UniMol array, inherited FDA query list `fda_drugs.csv`). The workflow never queries an external drug database and never generates drug features from SMILES.
- **Gene info** comes only from the XPert-author `pert_model/xpert/metadata/l1000_gene_info_978.csv` and the DToxS package. No external gene databases.
- The 1,004 screened drugs are a strict subset of the author's L1000 table; their embeddings are the author reference byte-identical (step 01 parity gate).

## Structure (finalized 11-step layout)

- `01_drug_embed/` — merged drug list (`drug_list.py`) + author-reference HG embeddings (`drug_embeddings.py`, byte-parity gate).
- `02_cm_baseline/` — consensus cardiomyocyte baseline (`cm_baseline.py`), authoritative gene mapping.
- `03_xpert_pretrain/` — the workflow-owned pretrain unit: used pretrain weight + `config/l1000.yaml` (data roots point at `pert_data/l1000/data/` and `input/`).
- `04_xpert_finetune_bench/` — XPert DTOXS-split benchmark (`xpert_finetune_bench.py`, `create_dtoxs_splits.py`).
- `05_xpert_finetune/` — all-train DTOXS finetune (`xpert_finetune.py`) + step-owned `config/dtoxs.yaml`.
- `06_xpert_predict/` — merged perturbed-transcriptome prediction (`perturb_transcriptome.py`) + z-score (`zscore.py`).
- `07_miri_phenotype/` — integrated MIRI reference gene sets + rank reference (`miri_phenotype.py`).
- `08_miri_efficacy_scores/` — the three efficacy scorers + z-components + plots (`miri_efficacy_scores.py`).
- `09_filter/` — WCS<0 + >=6/8 critical-pathway filters (`filter.py`).
- `10_integrate/` — integrated score over all drugs, ranked candidates, deterministic report (`integrate.py`).
- `11_report/` — agentic bilingual final drug report (`run_drug_agents.py`, `assemble_report.py`, `assemble_cn_report.py`, `prompt.md`; per-drug JSON cache in `cache/`, outputs in `output/`).
- `input/` — **real directory** owning the frozen UniMol array (`all_drugs_unimol_arr.npy`), idx→SMILES mapping (`all_drugs_idx2smi.npy`), DToxS h5ad, FDA drug list, curated `MIRI_pathways.json`, and the **author HG reference** (`HG_drug_embeddings.npy`).
- `_lib/utils.py` — workflow-root-aware paths; `N_L1000_DRUGS = 8981`; `resolve_l1000_drug_rows()`; `AUTHOR_HG_REFERENCE = input/HG_drug_embeddings.npy`.
- `_lib/scoring.py` — the three efficacy scorers + integration.
- `_lib/gene_mapping.py` — **authoritative DToxS→L1000 landmark gene map**, derived from **XPert-author + DToxS data only**. Resolves the 18 legacy-symbol landmarks (e.g. WASHC5 ← DToxS KIAA0196, ELP1 ← IKBKAP, JPT2 ← HN1L) using the author's own `l1000_gene_info_978.csv` `gene_name` column. `corrected_gene_name` is never a match key — it rejects the JPT2 ← CRAMP1 trap (CRAMP1 is Entrez 57585, a *different* gene; the DToxS table carries HN1L, JPT2's real alias, Entrez 90861).
- `_lib/verify_genes.py` — the gene-mapping gate (G1 completeness, G2 author-metadata alias verification, G3 frozen-h5ad audit), writes `_lib/gene_mapping_report.md`.

## Gene matching (why it matters)

DToxS (2015-era 3'-DGE symbols) contains all 978 L1000 landmark genes, but 18 landmarks only match via a **retired/alternative symbol**. Direct symbol matching zeroes them (v1's step-03 baseline lost 29 genes). The finalized baseline maps **978/978** landmarks and only the 11 genuinely undetected genes (SKIV2L, MIF, UBR7, PAK6, HSPA1A, BTK, PSMB8, PLEKHM1, HLA-DMA, PLSCR3, TSTA3) stay at 0.

## Key invariants (gates)

- **Step 01 parity** — `01_drug_embed/output/all_drugs_idx2HG.npy` is byte-identical to the author reference (`np.array_equal`).
- **Step 01 scope** — every drug in `01_drug_embed/output/drug_list.csv` has `drug_idx < 8981` (exact SMILES match). No FDA-extension rows are ever predicted.
- **Gene gates** — `_lib/verify_genes.py`: G1 all 978 landmarks resolve; G2 all 18 aliases author-metadata verified; G3 no confirmed h5ad gene misassignment. Any G1/G2 failure aborts (exit 1); G3 findings are reported.
- The DToxS finetune (steps 04/05) uses the shared frozen `input/dtoxs.h5ad` (not rebuilt).

## Entry points

```bash
cd workflow/xpert-fda-miri-v1
python 01_drug_embed/drug_list.py                # 1,004 L1000-contained FDA drugs
python 01_drug_embed/drug_embeddings.py          # author-reference HG (no HG training)
python _lib/verify_genes.py                      # gene-mapping gates (run before step 02)
python 02_cm_baseline/cm_baseline.py
# GPU steps (04 benchmark / 05 finetune / 06 predict) then 07–10
python 10_integrate/integrate.py
# agentic report (external biomni backend) — resumable
python 11_report/run_drug_agents.py
python 11_report/assemble_report.py
python 11_report/assemble_cn_report.py
python _lib/build_screen_notebook.py             # regenerate screen.ipynb
```

## Gotchas

- Never regenerate UniMol features from SMILES at predict time — step 06 indexes the precomputed array by `drug_idx`.
- The HG table is 8981×256 (author reference only); do not point the config at an 11104-row extension table.
- Pretrain weight + both YAMLs (pretrain + finetune) are workflow-owned; the model library resolves configs against the workflow root.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`), `11_report/cache/` JSONs, and `.env` are gitignored — never commit them.
- Steps 04/05 re-run the XPert DToxS finetune/benchmark on the author-reference HG — GPU-fp-nondeterministic across environments.
- The 11_report biomni agent backend (`/home/chang/cc-research`, overridable via `BIOMNI_CC_DIR`) is an external dependency; it is not vendored.
