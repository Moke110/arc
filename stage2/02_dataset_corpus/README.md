# 02_dataset_corpus — Stage 2 dataset corpus

**Status:** skeleton in build (ticket #14, 2026-11-04). The directory layout,
the code move and the Tahoe corpus move are done; this is **not yet the corpus
contract**. Tickets #15 (four CM corpora on native gene axes), #16
(seven crosswalk tables + coverage README) and #20 (docs and gates) fill the
rest. The governing decision is
[ADR 0033](../../docs/adr/0033-stage2-dataset-corpus.md).

## Purpose

One step owns every Stage 2 source dataset as an information-maximal h5ad
corpus: the substrate between a raw download and a model's training input. Model
adaptation (tiers, gene-space projection, author shims) belongs to the consumer
that builds a view, not here (ADR 0033 §5). `stage2/04_tahoe_data/` was folded
into this step and removed; `05`/`06` keep their numbers, so 04 is a gap in the
sequence, not a renumber.

From this decision the dataset ids are `dtoxs_primary`, `dtoxs_ipsc`, `sctox`
(was `screen_ctox`) and `tenv` (was `temp_env`). Frozen records keep their old
strings; the corpus README will carry the rename table (ticket #16).

## Layout

```
02_dataset_corpus/
├── AGENTS.md  README.md
├── execution/                        # code, intermediates, build records
│   ├── build_cm_datasets.py          # CLI today (ticket #15 replaces it)
│   ├── cmdata/                       # CM library; views.py retired
│   ├── tahoe/build_pseudobulk.py     # Tahoe corpus builder (moved from 04)
│   ├── _parts/                       # 14 per-plate CSR parts + _state (local)
│   ├── build_record/                 # tracked QA/manifests
│   └── tests/
└── output/                           # artifacts a later step reads
    └── corpus/v2/
        ├── tahoe.h5ad                # 67,018 × 62,710 CSR int32 (25.8 GB, local)
        └── tables/                   # seven crosswalk tables (ticket #16)
```

The four CM corpora (`dtoxs_primary.h5ad`, `dtoxs_ipsc.h5ad`, `sctox.h5ad`,
`tenv.h5ad`), the seven tables and `output/corpus/v2/README.md` do not exist
yet. No corpus is rebuilt from this skeleton state.

## What this step owns (target contract, ADR 0033)

`X` = summed raw counts (`int32`), `layers['log1p_cp10k'] = log1p(1e4·x/row_sum)`
(`float32`), native gene axis per source, a CM `obs` contract for the four CM
corpora and the Tahoe-native nine-column `obs` for `tahoe`, empty
`uns`/`obsm`/`varm`/`obsp`/`varp`, no integer encodings, no tier filter, no Δ.
The five corpora: `dtoxs_primary` 779 × 26,328 · `dtoxs_ipsc` 312 × 26,328 ·
`sctox` 3,072 × 31,880 · `tenv` 4,510 × 19,613 · `tahoe` 67,018 × 62,710.

## Moved in ticket #14

- `02_cm_data/execution/` → `execution/` (code; the five v1 views moved to
  `stage2/_frozen/cm-data-v1/` in ticket #13).
- `04_tahoe_data/execution/build_pseudobulk.py` → `execution/tahoe/`;
  `04_tahoe_data/execution/build_pretrain_view.py` →
  `05_prnet_pretrain/execution/` (the view is consumer-owned).
- `04_tahoe_data/output/pseudobulk/tahoe_pseudobulk.h5ad` →
  `output/corpus/v2/tahoe.h5ad` (renamed only; 25,805,562,874 B, shape
  67,018 × 62,710).
- `04_tahoe_data/output/pseudobulk/{manifest,qa,qa_report,qa_contrasts}` →
  `execution/build_record/tahoe_pseudobulk_*`.
- `04_tahoe_data/output/pseudobulk/{_parts,_state}` → `execution/_parts/`.
- `04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad` and its four
  tracked side files → `05_prnet_pretrain/output/pretrain_view/` (unchanged,
  `content_sha256 = 476812c0…`).
- Deleted: `04_tahoe_data/output/pretrain_view/_view_parts/` (8.26 GB gather
  cache; rebuild is 175 s) and the empty `02_cm_data/` + `04_tahoe_data/` trees.

Move evidence (hashes, bytes, inodes, disk before/after, `_view_parts`
reclaimed, `--check` peak RSS) is in
`execution/build_record/tahoe_corpus_move.json` and
`execution/build_record/tahoe_cm_schema_move.json`.

## Running today

From the repository root (the repo root is found by walking up to `.git`, never
by counting `../..` — ADR 0030 §4):

```bash
PY=/home/chang/miniconda3/envs/arc/bin/python
$PY stage2/02_dataset_corpus/execution/build_cm_datasets.py --check   # source counts
```

`--check` loads the four raw CM sources from `pert_data/` (measured 2026-11-04:
54.8 s wall, peak RSS 7,853,744 KiB ≈ 7.49 GiB); per the execution-environment
rule, run it under `/usr/bin/time -v` and do not run it alongside another heavy
job. The measured result is recorded in
`execution/build_record/tahoe_corpus_move.json`.
