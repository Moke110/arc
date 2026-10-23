# temp_env_all_5k — 5,000-gene HVG sensitivity benchmark

**Status (2026-10-23):** executed. 60/60 units completed (5 baselines + 7 real
models × 5 cold-drug folds). PRnet ranks first at **0.2735** macro AUPRC
(k=200), ahead of `delta_mean` at **0.2491**; it is the only entrant above the
baseline on this axis (+0.0243, paired drug CI +0.0113..+0.0367). This is a
sensitivity fork of [`../temp_env_all_18k/`](../temp_env_all_18k/), not a
replacement for it; the 18k archive and ADR 0027 are unchanged. Full analysis:
[`REPORT.md`](REPORT.md). Decision:
[ADR 0028](../../../docs/adr/0028-stage2-5k-hvg-benchmark.md).

## What this run changes

Only the gene axis. Wells, Δ, controls, groups, five cold-drug folds, entrants,
seeds, budgets, adapters and metrics are identical to the 18k run, and
`build_data.py` verifies that parity column-by-column before writing.

- **Axis:** the 5,000 genes with the smallest sum of per-source HVG ranks over
  `dtoxs_cm`, `screen_ctox` and `temp_env_all-smiles`.
- **Statistic:** Scanpy Seurat `dispersions_norm` on `layers["log1p_cp10k"]`,
  all wells, no `batch_key`.
- **Rank rule:** per source, finite mean in `[0.0125, 3]` → ranks `1..M` by
  descending dispersion (ties by corpus order); every other gene, including
  zero-mean/NaN, gets rank `18,165`. Sum the three ranks, keep the 5,000
  smallest (ties by corpus order), store in corpus order.
- **Prevalence:** k stays absolute at 100 / 200 / 500, so random prevalence is
  2.0% / 4.0% / 10.0% — five times the 18k denominators. The `no_change` row
  reports this floor (0.0400 at k=200).
- **Selection leakage:** the axis is unsupervised and uses all wells, including
  wells of drugs that are test in some fold. It is fixed across folds and
  declared in the ADR; no fold-specific gene axis exists.

Selection provenance is in `data/gene_set_5k.manifest.json`; the full per-gene
rank table is `data/gene_set_5k_ranks.csv`; the selected rows are
`data/gene_set_5k.csv`. The selection matches the pre-registered rule exactly
and is reproducible from the three source h5ads.

## Layout

| path | content |
|---|---|
| `data/` | frozen 5k artifacts: `counts`, `log1p_cp10k`, Δ, group Δ, controls, 5k true sets, identical splits, the gene-set tables, and `manifest.json` |
| `scripts/build_gene_set.py` | recomputes the 5k HVG selection and its provenance manifest |
| `scripts/build_data.py` | rebuilds the frozen data from the same h5ad, subsets the gene axis, recomputes true sets, and verifies 18k parity |
| `scripts/runner.py` | fold driver, subprocess bridge, checkpoint/resume, per-fold artifacts |
| `scripts/baselines.py` | train-mean, per-dose mean, ridge-on-fingerprint, nearest-drug baselines |
| `scripts/evaluate.py` | group/drug/fold/headline aggregation, bootstrap, `metrics.csv` / `per_drug.csv` / `leaderboard.csv` |
| `scripts/run_queue.sh` | 12 entrants × 5 folds under the shared GPU lock and 8 GB RAM scope |
| `scripts/tests/` | gene-set, data-contract, metric, harness, registry and adapter contract tests |
| `configs/models.json` | the same 12-entry registry as 18k |
| `models/<model>/` | author-default config, run card, and `FIDELITY.md`, re-declared for the 5k axis |
| `output/`, `results/` | per-unit artifacts and the evaluated leaderboard (gitignored) |

## Commands

```bash
# rebuild the selection and the frozen data (idempotent)
python scripts/build_gene_set.py
python scripts/build_data.py
python scripts/build_data.py --check

# offline contracts (no GPU)
make test

# five CPU baselines on fold 0
bash scripts/run.sh --dry-run

# one fold-0 smoke per real adapter
BENCH_SMOKE=1 bash scripts/run.sh --model prnet --fold 0

# the full 60-unit queue
bash scripts/run_queue.sh

# evaluate and summarise after the queue
python scripts/evaluate.py
python scripts/summarize_run.py
python scripts/aggregate_fidelity.py
python scripts/make_summary.py
```

Adapters run in the shared Stage 2 environments (`stage2/_envs/`), unchanged
from the 18k run.

## Fold-0 smoke results (5k)

Smoke runs use the adapters' shortened budgets (`BENCH_SMOKE=1`), one epoch or
a few steps. They prove wiring and shape, not ranking. All wrote finite
`(837, 5000)` predictions.

| model | fit wall (s) | total wall (s) | peak RSS (GB) | smoke AUPRC k=200 |
|---|---:|---:|---:|---:|
| PRnet | 5.0 | 16 | 1.1 | 0.1236 |
| CPA | 16.7 | 32 | 1.6 | 0.1228 |
| chemCPA | 9.4 | 18 | 1.8 | 0.1493 |
| biolord | 16.7 | 38 | 3.6 | 0.1240 |
| PerturbNet | 16.5 | 39 | 5.8 | 0.1548 |
| Squidiff | 10.0 | 21 | 1.6 | 0.0550 |
| CMonge | 71.2 | 89 | 1.4 | 0.1222 |

## Relationship to the 18k archive

`temp_env_all_18k/` remains the primary, full-transcriptome result. This run
answers one question: does replacing the 18,164-gene axis with a 5,000-gene
HVG axis change the ranking? It does: PRnet moves from −0.0345 to +0.0243
against `delta_mean`, while PerturbNet moves from −0.0128 to −0.0563. Absolute
AUPRC is not comparable across the two runs because prevalence rises from 1.1%
to 4.0%; use the paired gaps and ranks. Squidiff still fails to converge, so
the gene space was not the cause of its 18k failure.

## Results

- `results/leaderboard.csv` — 12 rows, k=200 headline, five folds.
- `results/metrics.csv`, `per_drug.csv`, `bootstrap.npz` — full evaluation.
- `results/comparison_18k.csv` — paired per-drug gap and its change vs 18k.
- `results/run_cards.csv`, `fidelity.csv` — execution and ledger aggregates.
- `summary.png`, `summary_table.csv` — metrics table and per-metric bars.
- `REPORT.md` — leaderboard, paired analysis, Squidiff failure and limitations.
