# 03_bench_models — Stage 2 benchmark harness

**Status:** harness, baselines, and all seven real adapters built (2026-10-21).
`data/` is the frozen ticket-11 build; `runner.py` + `evaluate.py` run the full
fit → predict → score path on CPU/GPU and write a leaderboard. All seven real
models (PRnet, CPA, chemCPA, biolord, PerturbNet, Squidiff, CMonge) completed a
one-fold smoke; ticket 16 runs the five cold-drug folds. This directory does not
import model internals — every adapter is a subprocess in its own environment.

## Quick start

```bash
# CPU smoke: five mandatory baselines, fold 0, leaderboard
bash stage2/03_bench_models/run.sh --dry-run

# one baseline / one or more folds
bash stage2/03_bench_models/run.sh --model ridge_fp --fold 0

# one real adapter, one fold (the smoke adds BENCH_SMOKE=1, see models/)
BENCH_SMOKE=1 bash stage2/03_bench_models/run.sh --model prnet --fold 0

# the ticket-16 queue: seven models × five folds under the GPU lock + RAM scope
bash stage2/03_bench_models/run_queue.sh
bash stage2/03_bench_models/run_queue.sh --force   # rerun completed units

# rerun completed units
bash stage2/03_bench_models/run.sh --dry-run --force

# offline contract + harness tests
cd stage2/03_bench_models && make test
```

`--dry-run` selects the five baselines from `configs/models.json`, fold 0, and
writes `output/dry_run/` + `results/dry_run/`. A normal run writes `output/`
and `results/`. `--model` and `--fold` are repeatable. Completed units are
skipped unless `--force` is passed.

## Layout

| path | content |
|---|---|
| `data/` | frozen ticket-11 artifacts (counts, Δ, groups, controls, true sets, splits); build with `build_data.py --check` |
| `build_data.py` | ticket-11 builder (not rerun by the harness) |
| `cases.py` | per-fold case bundles: train/valid with expression and Δ, test metadata only |
| `runner.py` | fold driver, subprocess bridge, checkpoint/resume, per-fold artifacts |
| `run_queue.sh` | ticket-16 seven-model × five-fold queue (GPU flock + 8 GB RAM scope) |
| `evaluate.py` | group/drug/fold/headline aggregation, bootstrap, `metrics.csv` / `per_drug.csv` / `leaderboard.csv` |
| `metrics.py` | metric primitives shared with the build sanity cases |
| `baselines.py` | train-mean, per-dose mean, ridge-on-fingerprint, nearest-drug baselines |
| `adapter_api.py` | `Adapter` base class and standard fit/predict CLI (carries `--controls`) |
| `adapters/no_change/` | CPU smoke adapter (`fit` no-op, `predict` zeros) |
| `adapters/baselines/` | CLI for the four non-trivial baselines |
| `adapters/prnet/`, `adapters/cpa/`, `adapters/chemcpa/` | ticket-13 real adapters |
| `adapters/biolord/`, `adapters/perturbnet/`, `adapters/squidiff/`, `adapters/cmonge/` | ticket-14/15 real adapters |
| `models/<model>/` | author-default config, run card, and `FIDELITY.md` per real model |
| `envs/` | `pip freeze` locks for the seven bench envs |
| `configs/models.json` | model registry (5 baselines + 7 real models) |
| `tests/` | ticket-11 contract tests, metric sanity, harness dry run, per-adapter offline contracts |
| `output/` | per-unit run dirs, shared `_cases/`, `run_status.csv` (gitignored) |
| `results/` | evaluation artifacts; `results/dry_run/` for the smoke run |
| `data/README.md` | data contract, shapes, hashes |
| `results/README.md` | metric artifact legend |

## Case bundles

The runner builds one shared bundle per fold under
`output/_cases/fold_<k>/`; all models read the same files.

- `train.npz` / `valid.npz`: one row per well; metadata plus `X` (int32),
  `log1p_cp10k` (float32), `delta` (float32), `well_row`, `delta_row`.
- `test.npz`: metadata only (`well_id`, `group_id`, `drug_key`,
  `compound_name`, `smiles`, `dose_um`, `dose_log`, `time_h`, `cell_type`,
  `plate`, `control_group`, `replicate`, `pharm_atc_moa`). It contains **no**
  `X`, `log1p_cp10k`, `delta`, or expression row index; `cases.verify_cases`
  fails closed on a violation.
- Controls are shared once via `data/controls.npz` and the case's
  `control_group`; the adapter receives `--controls` and the config's
  `data.controls` path.
- Test well IDs are visible (the adapter can know group sizes/plates); adapters
  must not join them back to the global expression arrays. The config no longer
  points at `group_delta.npy` or other truth arrays.

## Adapter seam

Each adapter is a subprocess with the ticket-05 CLI:

```bash
run.py fit     --train <npz> --valid <npz> --config <json> --run-dir <dir> [--controls <npz>]
run.py predict --cases <npz> --checkpoint <path> --out <npy> [--controls <npz>]
```

`fit` returns a checkpoint directory and writes `checkpoint_path.txt`; `predict`
returns `(n_cases, 18164)` float32. The registry env spec selects the
interpreter: `{"kind": "python"}` (current), `{"kind": "venv", "path": ...}`, or
`{"kind": "conda", "name": ...}`. No benchmark module imports model code.
Both subcommands also receive `--controls`; the real adapters use the shared
DMSO means as the control input and as the counterfactual basal state.

## Per-unit artifacts

`output/<model>/fold_<k>/` contains `config.json`, `run.log`,
`checkpoint/`, `pred_valid.npy`, `pred_test.npy`, `well_ids.json`,
`group_ids.json`, `group_pred_*.npy`, `group_true_*.npy`, `valid_metrics.json`,
`metrics.json`, and `status.json`. `output/run_status.csv` aggregates units.
A unit is resumed only when `status == completed` and all artifacts exist.
A timed-out fit predicts from its last checkpoint if present and stays
`timed_out`; a crash stays `failed` and never fabricates a prediction.

## Metrics

- AUPRC at k=100/200/500, direction macro, group macro; 200 is the headline.
- PCC, Spearman, R² (about the true mean), RMSE (raw `log1p_cp10k` Δ), NRMSE.
- `metrics.csv` uses `k=0` + `direction=macro` for regression rows;
  `auprc_up`/`auprc_down` are audit rows, never a second leaderboard.
- No micro AUPRC pooling (ticket 07 explicitly forbids it); `pooled_rmse` is
  the pooled regression view.
- `subset=atc-moa` is a descriptive view. A constant prediction has NaN
  PCC/Spearman (counted, labelled), not 0.
- Headline = mean of available fold macro AUPRCs at k=200. Five folds are
  required for `winner_eligible=true`; one-fold dry runs are `partial`.
- Bootstrap: 2,000 stratified 92-drug resamples, seed `20261022`, paired
  adjacent differences with percentile 95% CI.

## Ticket-13 real adapters (fold 0, 1-epoch smoke, 2026-10-21)

| model | env | fit / predict wall | test AUPRC k=200 | PCC | Spearman | R² | RMSE | checkpoint |
|---|---|---|---|---|---|---|---|---|
| `prnet` | `arc-stage2-bench-prnet` (py3.9, torch 2.0.1+cu118) | 4.3 s / 6.9 s | 0.0761 | 0.0206 | 0.1936 | −13.64 | 0.4266 | 4.7 M params |
| `cpa` | `arc-stage2-bench-cpa` (py3.9, torch 2.0.1+cu118) | 14.0 s / 10.8 s | 0.0759 | 0.0194 | 0.1488 | −15.14 | 0.4497 | 30.3 M params |
| `chemcpa` | `arc-stage2-bench-chemcpa` (py3.9, torch 2.8.0+cu128) | 7.7 s / 5.8 s | 0.0814 | 0.0665 | 0.2314 | −7.35 | 0.3392 | 14.6 M params |

These are one-epoch wiring proofs (`BENCH_SMOKE=1`), not the author recipes:
every fit stopped after the first epoch, and chemCPA is intentionally
un-pretrained. All three write finite `(837, 18164)` test predictions and 276
group means. Run cards, author defaults, and fidelity ledgers are in
`models/{prnet,cpa,chemcpa}/`.

## Ticket-14/15 real adapters (fold 0, tiny-budget smoke, 2026-10-21)

| model | env | smoke budget | fit / predict wall | test AUPRC k=200 | PCC | checkpoint |
|---|---|---|---|---|---|---|
| `biolord` | `envs/biolord` (py3.13, torch 2.12+cu130) | 1 epoch, 2 batches | 17 s / 16 s | 0.0764 | 0.0209 | 203 M params, 812 MB |
| `perturbnet` | `arc-stage2-bench-perturbnet` (py3.7, torch 1.13+cu117) | 1 epoch per stage | 21 s / 26 s | 0.0968 | 0.1570 | 169.8 M cINN, 743 MB |
| `squidiff` | `arc-stage2-bench-squidiff` (py3.13, torch 2.14+cu130) | 20 steps / 10 DDIM | 16 s / 9 s | 0.0224 | −0.00002 | 164.6 M, 628 MB |
| `cmonge` | `envs/cmonge` (py3.11, CPU JAX) | 4 AE steps + 20 CM iters | 81 s / 9 s | 0.0760 | 0.0199 | 19.2 M AE, 77 MB |

These are budget-truncated smoke runs, not the author schedules: biolord is at
1/500 epochs, PerturbNet at 1 epoch per stage, Squidiff at 20/100,000 steps,
CMonge at 20/10,000 iterations. Squidiff's 20-step sample is noise-scale by
construction (valid ε-MSE is still falling), so its near-zero PCC is expected.
The 276 group means are 92 test drugs × 3 doses. Run cards, author defaults,
and fidelity ledgers are in `models/{biolord,perturbnet,squidiff,cmonge}/`.

## Fidelity labels in the leaderboard

`ideal_fit` is now a string column: `full` for PRnet, CPA, Squidiff, and
CMonge; `partial` for chemCPA (978 pretrain omitted), biolord (transductive),
and PerturbNet (no control transcriptome); `none` for the CPU baselines. The
leaderboard has one ranking; the label is a report column, not a second table.

## Dry-run baseline scores (fold 0, 2026-10-21)

| model | headline k=200 | 95% CI | k=100 | k=500 | PCC | Spearman | R² | RMSE |
|---|---|---|---|---|---|---|---|---|
| `delta_mean` | 0.214468 | 0.176494–0.253793 | 0.199932 | 0.235665 | 0.2536 | 0.2073 | 0.0295 | 0.1403 |
| `delta_mean_dose` | 0.209090 | 0.171911–0.247085 | 0.194796 | 0.230006 | 0.2508 | 0.2025 | 0.0256 | 0.1404 |
| `ridge_fp` | 0.101021 | 0.080352–0.122417 | 0.088673 | 0.125705 | 0.0871 | 0.0590 | −0.7501 | 0.1718 |
| `nearest_drug` | 0.081212 | 0.065257–0.098796 | 0.068459 | 0.105523 | 0.0708 | 0.0506 | −1.7145 | 0.1992 |
| `no_change` | 0.011011 | 0.011011–0.011011 | 0.005505 | 0.027527 | NaN | NaN | −0.0210 | 0.1450 |

`no_change` reproduces the prevalence exactly at k=200; its NaN PCC/Spearman is
the contract's constant-prediction policy. These are one-fold smoke values, not
the final five-fold ranking.

## Next

Tickets 14–15 are complete: all seven real adapters are registered, smoke-tested
on fold 0, and carry a config, run card, and fidelity ledger. Ticket 16 runs the
five-fold queue for the seven entrants plus the five CPU baselines.
