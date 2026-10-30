# Tahoe-100M data home and the Stage 2 bulk pretrain corpus

Two things live here: the **whole atlas** under `pert_data/Tahoe-100M/source/` and the
**well-level pseudobulk corpus** built from it, per ADR 0029.

## The atlas (`pert_data/Tahoe-100M/source/`)

Imported on 2026-10-30 from a Windows-side HuggingFace snapshot
(`E:\Desktop\tahoe\data\hf`) rather than re-downloaded: 3,388 expression shards
(95,624,334 cells), the five metadata tables, and the authors'
`pseudobulk_differential_expression` (1,026 shards, 4,089,820,780 rows) —
428.8 GB total. Provenance and the verification record are in
`pert_data/Tahoe-100M/local/Tahoe_summary.md`.

```bash
# offline verification (no network needed)
python pert_data/Tahoe-100M/script/verify_import.py

# authoritative sha256 against HF's published hashes (needs network)
python stage2/04_tahoe_pretrain/script/fetch_tahoe.py --verify-only deep \
    --scope full --include-pseudobulk
```

`fetch_tahoe.py` remains the fetch/verify tool; `--scope full` is an in-place
top-up, and the 36-shard vendored subset it originally produced is a strict
subset of what is on disk now.

Final state (2026-10-30): the authoritative deep check **passed** — 4,425 files
compared against HF's published `lfs.sha256`, `ok: true`, 0 problems
(`logs/tahoe_verify_deep_1841.log`). `MANIFEST.jsonl` now holds those 4,425 HF
hashes, so the offline check above is anchored to HF. Provenance and reading
notes: `pert_data/Tahoe-100M/local/Tahoe_summary.md`.

## The corpus (`output/pseudobulk/`)

Built by `script/build_pseudobulk.py` in resumable stages:

```bash
python script/build_pseudobulk.py scan        # shard -> plate(s) map  (~20 s)
python script/build_pseudobulk.py aggregate   # 14 per-plate CSR parts (~50 min)
python script/build_pseudobulk.py assemble    # sparse h5ad, streamed  (~4 min)
python script/build_pseudobulk.py qa          # vs the authors' DESeq2 DE
```

| | |
|---|---|
| unit | `(plate, well, cell_line)` — 67,018 profiles |
| gene space | Tahoe native 62,710 (`token_id - 3`) |
| `X` | summed raw UMI counts (int32) |
| `layers['log1p_cp10k']` | `log1p(1e4 * x / row_sum)` |
| controls | 1,400 `DMSO_TF` profiles, `control_group = tahoe::plateN` |
| cells consumed | 95,624,334 — **every cell exactly once** (matches the card) |
| file | `tahoe_pseudobulk.h5ad`, 25.8 GB |

Δ is deliberately **not** materialized: it is
`log1p_cp10k(treated) − mean(same-plate, same-cell-line DMSO_TF)`, one groupby
away, and freezing it would freeze a pairing rule that belongs to the model.

### Two traps this build encodes

1. **A well is not a bulk profile.** Every well pools most of the 50 cell lines,
   so summing a well gives a mixture whose DE is a Simpson's-paradox artefact.
   The corpus is cell-line-stratified; `(plate, well, cell_line)` is the unit.
2. **A shard is not always one plate.** 13 of the 3,388 shards straddle a plate
   boundary (e.g. shard 236 holds plate4 + plate5 cells), so aggregation routes
   each *cell* by its own `plate` column and reads the shard once per plate it
   touches.

### QA

`output/pseudobulk/qa_report.md` + `qa.json` + `qa_contrasts.csv`. Gates:
median detected genes ≥ 12,000 per profile; every `(plate, cell_line)` has a
vehicle; median per-contrast Pearson ≥ 0.60 against the authors' DESeq2
`log2FoldChange` on genes with `baseMean > 500`. The ADR explains why the
originally proposed "Spearman ≥ 0.8 over all genes" was retired after
measurement: over the full gene set the statistic is ~0.2 because near-null
genes dominate ranks, and even on well-measured genes the achievable agreement
is capped by the fact that the shipped DE was produced under a
normalization/aggregation that is not a naive count sum.

Partial rebuilds: `--plates plate4,plate5` limits any stage, and each plate's
part is skipped if it already exists (`--force` to rebuild).
