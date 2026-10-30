# 04_tahoe_data — Tahoe-100M pseudobulk pretrain corpus

The **whole-transcriptome pretraining input** for Stage 2: a bulk, well-level
corpus derived from Tahoe-100M and built once, then read by the pretrain /
adaptation work. CM data is bulk (one physical library per well), so the
pretraining substrate is pseudobulk profiles rather than single cells;
[ADR 0029](../../docs/adr/0029-tahoe-pseudobulk-pretrain-corpus.md) freezes how
they are built. The corpus is model-agnostic: identity (genes, units) and raw
counts are frozen; thresholds, Δ pairing and gene-space projections belong to a
view or to the model.

**This step performs no download.** The atlas is an external item
([`pert_data/Tahoe-100M/`](../../pert_data/Tahoe-100M/) and ADR 0030); its
importer and verifiers live in `pert_data/Tahoe-100M/script/` (`fetch_tahoe.py`,
`tahoe-dl.sh`, `run_fetch.sh`, `verify_import.py`, offline `verify.sh`). This
step only aggregates bytes that are already on disk.

## Layout

```
04_tahoe_data/
├── AGENTS.md  README.md
├── execution/build_pseudobulk.py      the four resumable stages
└── output/pseudobulk/
    ├── tahoe_pseudobulk.h5ad          25.8 GB corpus (local-only)
    ├── manifest.csv                   67,018 profiles x plate/well/line/drug/dose/n_cells
    ├── qa.json  qa_report.md          the gates and their measured values
    ├── qa_contrasts.csv               3,755 per-contrast agreement rows
    ├── _parts/                        14 per-plate CSR parts (resumable intermediate)
    └── _state/shard_plate.json        scan bookkeeping (machine state)
```

## Corpus

| | |
|---|---|
| unit | `(plate, well, cell_line)` — well is the upstream `sample`; **67,018 profiles** |
| gene space | Tahoe native **62,710** (`token_id - 3`), not the 18,164 CM intersection |
| `X` | summed raw UMI counts (`int32`) |
| `layers['log1p_cp10k']` | `log1p(1e4 * x / row_sum)` |
| controls | 1,400 `DMSO_TF` profiles, `control_group = tahoe::plate` |
| cells consumed | 95,624,334 — every cell exactly once (matches the atlas card) |
| file | `tahoe_pseudobulk.h5ad`, 25.8 GB |

Δ is deliberately **not** materialized: it is
`log1p_cp10k(treated) − mean(same-plate, same-cell-line DMSO_TF)`, one groupby
away, and freezing it would freeze a pairing rule that belongs to the model.

## Reproduce

Step-local program: `python`, with `numpy`, `pyarrow`, `scipy`, `h5py`,
`anndata`, `pandas` available (the build ran on the shared `arc` interpreter).
Run from the repository root; the script resolves the root by walking up to
`.git`, so any working directory works.

```bash
# report the resolved root / Tahoe input / output dir, then exit
python stage2/04_tahoe_data/execution/build_pseudobulk.py --check

# the four stages (resumable, in order)
python stage2/04_tahoe_data/execution/build_pseudobulk.py scan        # shard -> plate map  (~20 s)
python stage2/04_tahoe_data/execution/build_pseudobulk.py aggregate   # 14 per-plate parts  (~50 min)
python stage2/04_tahoe_data/execution/build_pseudobulk.py assemble    # sparse h5ad + manifest  (~4 min)
python stage2/04_tahoe_data/execution/build_pseudobulk.py qa          # vs the authors' DESeq2 DE
```

`--plates plate4,plate5` limits any stage; an existing plate part is skipped
(`--force` rebuilds it); `--shards N` limits scan/aggregate to the first N
shards for a smoke run. `qa` exits non-zero if any gate fails, so it doubles as
the corpus acceptance test.

## QA (2026-10-30 build)

From [`output/pseudobulk/qa.json`](output/pseudobulk/qa.json) and
[`qa_report.md`](output/pseudobulk/qa_report.md):

- **67,018 profiles**, of which **61,767 (92.2%) have ≥ 100 cells**.
- median detected genes per profile **25,290** — gate ≥ 12,000: **PASS**.
- **3,755 contrasts** compared against the authors' shipped DESeq2
  `log2FoldChange`:
  - `baseMean > 800` (median 328 genes/contrast): Pearson **0.812** (p10 0.669) — reported, not gated;
  - `baseMean > 500` (644 genes/contrast): Pearson **0.779**, Spearman **0.713** — gate ≥ 0.60: **PASS**;
  - `baseMean > 300`: Pearson **0.741**; all-gene Spearman **0.173** — reported only.
- **700 `(plate, cell_line)` pairs** carry a vehicle control; none missing — **PASS**.

The all-gene Spearman is context, not a gate: on genes with no real effect both
estimates are noise, so rank agreement there is uninformative. The originally
proposed "Spearman ≥ 0.8 over all genes" gate was retired after measurement
(full reasoning in ADR 0029 §QA): even on well-measured genes the achievable
agreement is capped because the shipped DE was produced under a
normalization/aggregation that is not a naive count sum.

## What is tracked, what is local

| path | git |
|---|---|
| `execution/build_pseudobulk.py`, `AGENTS.md`, `README.md` | tracked |
| `output/pseudobulk/{qa.json,qa_report.md,qa_contrasts.csv,manifest.csv}` | tracked build record (force-added past the `output/` and `*.csv` rules) |
| `output/pseudobulk/tahoe_pseudobulk.h5ad` | local-only, 25.8 GB, ignored |
| `output/pseudobulk/_parts/` | local-only resumable cache, ignored |
| `output/pseudobulk/_state/` | local-only machine state, ignored |

`_parts/` and `_state/` are **intermediates, not results** — they exist only so
that a failed or partial rebuild can resume. `_state/shard_plate.json` is the
scan's map from expression shard to plate(s); deleting it is safe, `scan`
recreates it. Both trees are hidden from git by the repository-wide `output/`
rule, which is why they need declaring here:

```console
$ git check-ignore -v stage2/04_tahoe_data/output/pseudobulk/_state/shard_plate.json
.gitignore:62:output/	stage2/04_tahoe_data/output/pseudobulk/_state/shard_plate.json
```

## Two traps this build encodes

1. **A well is not a bulk profile.** Every well pools most of the 50 cell lines,
   so summing a well gives a mixture whose DE is a Simpson's-paradox artefact.
   The corpus is cell-line-stratified; `(plate, well, cell_line)` is the unit.
2. **A shard is not always one plate.** 13 of the 3,388 shards straddle a plate
   boundary (e.g. shard 236 holds plate4 + plate5 cells), so aggregation routes
   each *cell* by its own `plate` column and reads the shard once per plate it
   touches.

## Scope

Pseudobulk fixes sparsity, not biology: this is a cancer-cell-line corpus at a
single 24 h timepoint with a 24/525 drug-key overlap against the CM corpus
(ADR 0029 §Consequences). It is a pretraining substrate, not a Stage 2
evaluation dataset. Atlas provenance, the HF revision and the verification
record are in
[`pert_data/Tahoe-100M/local/Tahoe_summary.md`](../../pert_data/Tahoe-100M/local/Tahoe_summary.md).
