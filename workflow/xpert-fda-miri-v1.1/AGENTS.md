# AGENTS.md

The **DToxS-LINCS testbed MIRI screening workflow** (`workflow/xpert-fda-miri-v1.1/`;
see docs/adr/0012): a 12-step mirror of the finalized
`workflow/xpert-fda-miri-v1/` (ADR 0007) whose only substantive change is the
XPert finetune data. Instead of v1's frozen `input/dtoxs.h5ad` (521 samples /
36 drugs, built outside the repo from the former `pert_data/dtoxs-xpert`
well-level deposit), steps 05/06/07 train on the **workflow-rebuilt
DToxS-LINCS h5ad** (`03_cm_finetune_data/output/dtoxs_lincs.h5ad`, 662 samples
/ 45 drugs) built from the **LINCS LDS-1593 well-level package** released flat
at `pert_data/dtoxs_primary/Data/Level3` on 2026-08-18 (the GSE146095/GSE146096
raw GEO deposit was deleted the same day as a strict subset; the v1-era
package is frozen byte-identically at `pert_data/dtoxs_primary/v1-frozen/`).
Everything else — the screened drug scope (1,004 L1000-contained FDA drugs,
`drug_idx < 8981`), step-02 baseline, MIRI phenotype/scoring, filters,
integration, and the agentic report — is v1-identical, preserved by the
12-step renumbering (same layout as xpert-fda-miri-v2). **v1 is untouched and
frozen.**

## Why v1.1 exists (testbed)

v1.1 tests two things before anything is ever changed in v1:
1. the **fresh DToxS-LINCS build** (step 03) and
2. the **step renumbering** to the 12-step layout.

It is not (yet) a canonical screen: its benchmark/report numbers differ from
v1's because the rebuilt finetune data is a legitimate superset (the LINCS
package carries the non-KI perturbations and extra anthracyclines/cardiac
drugs the raw GEO deposit lacks). v1.1's outputs are for comparison, not for
replacing v1.

## Scope and provenance (author + LINCS LDS-1593 data only)

- Screens **only** the FDA drugs already present in the L1000 embedding table
  (exact canonical-SMILES match, `drug_idx < 8981`; 1,004 of 3,127 drugs).
  Because every screened drug uses the author-reference HG embedding directly,
  **no HG training / low-LR finetune / FDA-extension inference is ever run**.
- **Drug info** comes only from the workflow-owned `input/` files (author L1000
  SMILES/HG tables, precomputed UniMol array, inherited FDA query list) and
  the LINCS small-molecule metadata (`pert_data/dtoxs_primary/Metadata/
  Small_Molecule_Metadata.txt`). The workflow never queries an external drug
  database and never regenerates drug features from SMILES.
- **Gene info** comes only from the XPert-author
  `pert_model/xpert/metadata/l1000_gene_info_978.csv` and the DToxS-primary
  gene index. No external gene databases.
- `pert_data/` is **read-only**: step 03 never modifies a source file; the
  DToxS-LINCS h5ad + maps + manifest + parity delta are regenerated workflow
  artifacts.

## Structure (12-step layout, mirrored from v2)

- `01_drug_embed/` — merged drug list (`drug_list.py`, 1,004 drugs) + author-reference HG embeddings (`drug_embeddings.py`, byte-parity gate). Identical to v1.
- `02_cm_baseline/` — consensus baseline from the DToxS-primary control wells of the raw GEO deposit (`pert_data/dtoxs_primary`, GSE146096); byte-identical to v1's baseline.
- `03_cm_finetune_data/` — **the v1.1 change**: rebuilds the DToxS finetune h5ad `output/dtoxs_lincs.h5ad` from the LINCS LDS-1593 well-level package (read-only) + the compound→drug-index manifest (`output/drug_manifest.csv`), the gene map (`output/gene_maps/`), the QN reference, the parity delta vs v1's frozen `input/dtoxs.h5ad` (`output/parity_delta.md`), and the composition report `output/finetune_report.html`.
- `04_xpert_pretrain/` — the workflow-owned pretrain unit (weight + `config/l1000.yaml`), carried over untouched from v1.
- `05_xpert_finetune_bench/` — DToxS benchmark on the **unified folds** (ADR 0014): `create_dtoxs_splits.py` writes `warm-start_fold-1..10` (drug-overlap) + `cold-drug_fold-1..10` (drug-disjoint) + `dtoxs_full` on the rebuilt h5ad (invariant gate, seed 2024); `xpert_finetune_bench.py` benches the 20 folds and overlays per-fold dots, and — as the **first column** (`pretrain`, amber) — zero-shot-evaluates the step-04 L1000 pretrain weight on all 20 fold test sets (no finetuning) under the workflow-owned `config/dtoxs_pretrain_eval.yaml` (pretrain architecture/bin `[0, 15.002]`; the drifted aux `class_fc` head is dropped by the library's tolerant test-mode load).
- `06_xpert_finetune/` — all-train DToxS-LINCS finetune (`xpert_finetune.py`, `config/dtoxs.yaml`) -> the final prediction model `output/dtoxs_dtoxs/`.
- `07_xpert_predict/` — perturbed-transcriptome prediction + z-score (`perturb_transcriptome.py`, `zscore.py`).
- `08_miri_phenotype/` — integrated MIRI reference gene sets + rank reference.
- `09_miri_efficacy_scores/` — the three efficacy scorers + z-components + plots.
- `10_filter/` — WCS<0 + >=6/8 critical-pathway filters.
- `11_integrate/` — integrated score over all drugs, ranked candidates, deterministic report.
- `12_report/` — agentic bilingual final drug report (external biomni backend).
- `input/` — frozen assets (author HG reference, UniMol array, idx→SMILES, FDA drug list, curated `MIRI_pathways.json`), hard-linked from v1.
- `_lib/utils.py` — workflow-root-aware paths; `CM_DATA_ROOT = pert_data`, `CM_DATASETS = ["dtoxs_primary"]`, `DTOXS_LINCS_ROOT`, `N_L1000_DRUGS = 8981`, `V1_DTOXS_H5AD` (parity reference).
- `_lib/gene_mapping.py` — DToxS-primary → L1000 landmark map (symbol space, legacy-symbol aware) + `load_dtoxs_lincs_raw()` (LINCS well-level reader).
- `_lib/drug_resolution.py` — compound→L1000 drug index (canonical-SMILES match, else author `cmap_name` name match) + the auditable manifest.
- `_lib/verify_genes.py` — the DToxS gene gates (G1 completeness, G2 alias verification, G3 drug-level audit), writes `_lib/gene_mapping_report.md`.

## The DToxS-LINCS finetune data (step 03)

- Reads the **LINCS LDS-1593 well-level package** (`pert_data/dtoxs_primary/
  Data/Level3`, the DToxS center's Level-3 deposit, 956 well records across
  Sequencing Experiments 1/2/4/5{Cor,Spec-Fixed,Spec-Var-N}),
  not the raw GEO deposit (which was deleted 2026-08-18 and held only the
  23 KI + CTRL DGE wells).
  Compound codes resolve via the LINCS `Small_Molecule_Metadata.txt`
  (`sm_center_batch_id` → name + SMILES).
- Only samples whose compound resolves to an HG drug (`drug_idx ∈ [0, 8981)`)
  enter; antibodies (TRS), proteins (IGF/TNF), and unmatched codes (END, DIC,
  SUN+) are excluded and recorded in the manifest. The rebuilt set is **662
  samples / 45 drugs** (vs v1 frozen 521 / 36) — a superset (non-KI drugs:
  estradiol, amiodarone, cyclosporin-A, milrinone, ...; anthracyclines;
  cardiac drugs) with 5 compounds re-indexed by the author canonical-SMILES
  rule (cyclosporin 73→72, crizotinib 4786→1176, estradiol 5674→1291,
  paroxetine 6743→6737, ruxolitinib 6798→6797).
- The h5ad is schema-identical to v1's `input/dtoxs.h5ad` in every field the
  XPert loader consumes (obs `pert_id, pert_idx, pert_dose[_idx], pert_time[_idx],
  cell_iname, cell_idx, tissue_idx` + additive `dataset_id`; var = 978 landmarks
  in L1000 `gene_idx` order; `obsm["X_ctl"]` per-sample matched control; X float32
  log2(CPM+1)). The DToxS single-dose experiment trains at `pert_dose 0 →
  dose_idx 0`, 24 h slot (`time_idx 2`), exactly like v1.
- Per-group **quantile alignment** to the L1000 pretrain pooled distribution
  (before matched-control computation); the finetune config bin range is
  `[0, 15.002]` (the pretrain range, ADR 0010).
- **Parity gate (ADR 0012)** — `03_cm_finetune_data/build_dtoxs_lincs.py` step
  4: requires **schema parity** (structure-identical to v1's `input/dtoxs.h5ad`
  in the loader fields) and writes an **auditable delta** report
  (`output/parity_delta.md`) of the sample/drug differences. Sample-for-sample
  identity with v1's frozen data is neither required nor enforced.

## Finetune-data composition report (step 03)

`03_cm_finetune_data/finetune_report.py` emits a self-contained HTML report
(`output/finetune_report.html`):
- a **collapsible drug list** grouped by the author L1000 `moa` (from
  `pert_model/xpert/metadata/l1000_sdst_drugs_info_8432.csv`), with
  `Unclassified` for compounds without a recorded moa;
- 2D **PCA + UMAP of cell embeddings** — L1000 pretrain aggregated per cell
  line (mean embedding) + every finetune sample, colored by cell type with the
  CM samples split by `dataset_id`. The embedding is the pretrained XPert
  `cell_Embeddings` output: binned expression (128 bins over [0, 15.002]) +
  per-gene PPI vectors, LayerNorm, mean-pooled over the 978 landmarks.

## Key invariants (gates)

- **Step 01 parity** — `01_drug_embed/output/all_drugs_idx2HG.npy` is
  byte-identical to the author reference (`np.array_equal`); `drug_list.csv` is
  1,004 rows (identical to v1).
- **Step 02 parity** — `02_cm_baseline/output/cm_ctl_landmark.csv` is
  byte-identical to v1's baseline.
- **Schema parity** — step 03: the rebuilt h5ad is structure-identical to
  v1's `input/dtoxs.h5ad` in the loader fields (asserted by the parity gate).
- **Gene gates** — `_lib/verify_genes.py`: G1 all 978 landmarks resolve, G2 all
  18 legacy aliases author-metadata verified, G3 drug-level audit. Any
  G1/G2 failure aborts (exit 1).
- **Drug-inclusion gate** — every obs `pert_idx` lies in `[0, 8981)`; the
  manifest classifies every compound (in_hg / extension / unmatched) and the
  exclusions are auditable.

## Entry points

```bash
cd workflow/xpert-fda-miri-v1.1
python 01_drug_embed/drug_list.py                # 1,004 L1000-contained FDA drugs
python 01_drug_embed/drug_embeddings.py          # author-reference HG (no HG training)
python _lib/verify_genes.py                      # gene gates (G1/G2 anytime; G3 needs step-03 h5ad)
python 02_cm_baseline/cm_baseline.py             # DToxS-primary consensus baseline
python 03_cm_finetune_data/build_dtoxs_lincs.py  # DToxS-LINCS h5ad + manifest + parity delta
python 03_cm_finetune_data/finetune_report.py    # composition report (HTML)
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
- The finetune/bench config (`06_xpert_finetune/config/dtoxs.yaml`) uses
  dataset name `dtoxs` (the library's `dtoxs` branch) with
  `DTOX_data_root` pointing at `03_cm_finetune_data/output/dtoxs_lincs.h5ad`.
- `create_dtoxs_splits.py` writes the unified fold columns — `warm-start`
  (10, drug-overlap) + `cold-drug` (10, drug-disjoint) + `dtoxs_full`
  (all-train) — on the rebuilt h5ad, with a built-in invariant gate (drug
  overlap 1.0/0.0, test-fraction band, coverage, no stale columns; ADR 0014).
- Step 07 predicts at `--time-idx 2` (24 h under the author time map).
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`), `12_report/cache/` JSONs,
  and `.env` are gitignored — never commit them.
- Steps 05/06 re-run the XPert finetune on the DToxS-LINCS data —
  GPU-fp-nondeterministic across environments. The `pretrain` benchmark column
  re-evaluates the step-04 L1000 weight (no finetuning) and is deterministic
  per GPU; the DToxS-specific embedding tables (dose/time/position) it drops
  are untrained in that column, exactly as they are at finetune epoch 0.
- The 12_report biomni agent backend (`/home/chang/cc-research`, overridable
  via `BIOMNI_CC_DIR`) is an external dependency; it is not vendored.
- v1 stays frozen: nothing in this workflow writes into v1 or modifies v1's
  inputs. The `input/` assets are hard-linked from v1 (read-only).
