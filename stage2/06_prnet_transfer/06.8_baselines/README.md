# 06.8 · baselines — five-fold simple baselines (ticket #8)

Status: **code finalized, fixture-verified, and the real five-fold record produced
2026-11-02** on the frozen 06.3 input (integration `e0e2d68`, pinned
`arc-stage2-bench-prnet` run; §3 measures it, §5 records it). No science result exists here:
this step produces the paired-comparison reference for #9, it does not show that any model
beats it.

Frozen input: ticket #3's `06.3_input/output/` — 118 main-set drugs / 951 wells, five
drug-level folds (per fold ≈ 570 train / 190 valid / 190 test wells; 16 plates). The truth Δ
and the direction AUPRC@200 come from `06.4_pipeline` only.

## 1. What it computes

| baseline | role | rule (fit uses **only** the fold's `role=train` treated responses) |
|---|---|---|
| `no_change` | **primary** | Δ = 0 for every well, i.e. predicted expression = the fold's fixed `c*` |
| `delta_mean` | **primary** | within-condition replicate mean first (condition = `drug_key\|dose_µM`), then equal-weight mean over the fold's training conditions; a scored condition whose exact key is in training uses that training mean, otherwise the pooled training mean (documented fallback; the run counts collisions) |
| `delta_mean_dose` | auxiliary | training condition means pooled per dose (`:g` µM); unseen dose → pooled training mean |
| `ridge_fp` | auxiliary | L2 ridge (`alpha=1.0`, fixed a priori, never tuned; `fit_intercept=True`, `solver="cholesky"`) from the study's own dose encoding (1024-bit rFCFP4 × `log10(dose_µM+1)`, the model's `drug_dim` label) to the training condition-mean Δ |

Predictions are written in the `ptlib.predict` format (`pred_expr.npy`, `delta.npy`,
`meta.csv` with the row ids, `manifest.json`, `Δ = prediction − c*`) and scored by
`ptlib.evaluate.run_evaluation` — **the one evaluator of the study**, never re-implemented.
The auxiliary baselines are listed separately in the manifest/config and never enter the
primary `no_change`/`delta_mean` comparison.

Leakage rule: only the fold's **training** treated responses enter a fit. Validation/test
treated rows are never read by a fit; the matched DMSO means used to form the training Δ are
controls and are never split by drug. `execution/tests/test_leakage.py` poisons every
held-out (valid / test) target in a copy of the fixture and requires fitted parameters and
predictions to stay byte-identical.

## 2. Record (`output/`, tracked)

| file | content |
|---|---|
| `summary_table.csv` | one row per fold × baseline plus a macro-over-folds row: headline AUPRC@200, up/down, RMSE, Pearson, Spearman |
| `metrics.csv` | the evaluator's long metric rows, keyed by fold × baseline |
| `per_drug_dose.csv` | per-fold, per-baseline, per-condition predictions/metrics (`condition`, `drug_key`, `dose_um`, `n_replicates`, direction set sizes, auxiliary metrics) — #9's paired-comparison input |
| `fold_coverage.csv` | fold sizes, row/`c*`/gene-axis hashes, `test_rows_disjoint_from_train`, conditions scored exactly once |
| `manifest.json` | schema `stage2/06.8_baselines/v1`, ticket, command, git SHA/branch, frozen-input artifact sha256 + manifest sha256 + locked-config sha256/mapping, gene axis, config (baselines, ridge alpha, folds, metric, pooling, memory policy), fold coverage, invariants, per-file record sha256, environment versions, timings, **peak RSS** |

Every fold also writes `folds/fold<k>/fold_record.json`, `<baseline>/{pred,eval}/` (the
pipeline-format prediction directory + the evaluator record) and `fits/<baseline>/`
(`fit_record.json` with the train-row/content hashes, `params.npz`). Those are local-only;
only the five files above are tracked (narrow `.gitignore` exception).

## 3. Memory policy and measured peaks (ticket #8)

The real-data row path was genuinely unbounded in the crashed-session version: it read the
whole `n_train_rows × 18,164` log1p layer, built a same-size matched-control matrix and a
float64 row-Δ matrix, and hashed a float32 copy of it. The rewrite removes every
`observation × 18,164` array:

- the dense CM log1p layer is read and aggregated **one `GENE_BLOCK` (1024) of genes at a
  time** (`iter_layer_blocks`); the largest temporary is `n_train_rows × 1024` float64;
- matched controls are kept as `n_plates × n_genes` plate means (`plate_control_matrix`),
  byte-identical to `CMControlPool.mean_for_rows`, and indexed per gene block;
- condition means are accumulated in float64 and stored as `n_conditions × n_genes`
  **float32**; `train_delta_sha256` hashes that condition-level Δ, so it is block-size
  independent;
- predictions are condition-level float32; the two `.npy` files are created with
  `np.lib.format.open_memmap` and filled in `ROW_BLOCK` (64) row blocks — no row-level
  prediction array in memory;
- one process, no worker pool (`n_jobs` is never set).

Measured (WSL2, `arc` env, `/usr/bin/time -v` inside `systemd-run --user --scope`; the
`06.4` fixture is 204 observations / 18,164 genes, the row-replicated fixture repeats every
row 18× = 3,672 observations = **3.4× the real per-fold train rows and 3.4× the real test
rows**, so it is a conservative upper bound):

| job | pre-fix | fixed |
|---|---|---|
| fit-table path alone (`build_fit_table`, 5 folds, 1,944 train rows) | 1,843 MB peak (983.9 MB after imports, **+859 MB**) | 1,084 MB peak (**+101 MB**) |
| `make dry-run` on the standard fixture, `MemoryMax=4G` | 1,261,316 kB (1.203 GiB) | 1,244,788 kB (1.187 GiB) |
| full `--all` on the 3,672-row replicated fixture (648 test rows/fold) | 1,922,484 kB (1.833 GiB) | 1,807,400 kB (1.724 GiB) |
| pinned-env (`arc-stage2-bench-prnet`) fold-0 smoke on the fixture | — | 998 MB |
| `make check` on the real 06.3 input (confirms the #4 header-read loader fix) | — | 881 MiB, 21.3 s |
| **real `make run` on 06.3** (5 folds, 951 test rows, 16.4 GB input hash) | — | **1,081 MB peak** (`/usr/bin/time -v` 1,107,292 kB = 1.056 GiB), 5 m 37 s |

The fixed full-run peak is now dominated by imports (~1 GB torch + numpy) and by the shared
evaluator's own dense test-row arrays (`pred_expr`, `delta`, float64 truth/pred/Δ,
`aggregate_rows` copies ≈ 6–7 × `n_test_rows × 18,164 × 4 B`), not by the baselines; at the
real ~190 test rows/fold that term is small. `make mem-check` reproduces the row-replicated
number (`systemd-run … -p MemoryMax=4G -p MemorySwapMax=4G -- make mem-check`).

**Correction to the crash record.** The 9.7 GiB / `4min 33s` scope in
[`docs/agents/execution-environment.md`](../../../docs/agents/execution-environment.md) was
**not** a baselines dry run. Journal boot `-2`: the scope started `2026-11-02T12:18:13` as
`bash -c "cd /home/chang/arc && (script/verify_all.sh > /tmp/verify_all_after_12.log …)"`
and ended `12:23:09` with `Consumed 4min 33.501s CPU time, 9.7G memory peak` (52 GB of
external-data hashing, mostly page cache). `baselines.py --all --data-dir
/tmp/t08-baselines-fixture` started `12:33:34`, and the VM died at `12:33:50` before the
scope logged a memory peak (`23.160s CPU time` only). The memory hardening in this step is
still required by the memory rule and is measured above.

## 4. Commands

```bash
cd stage2/06_prnet_transfer/06.8_baselines/execution
make test          # offline unittest suite (arc env, fixture): 38 tests, 5 files
make check         # contract check against the frozen 06.3 input, writes nothing
make dry-run       # full five-fold record on the synthetic fixture (pipeline check, not evidence)
make mem-check     # row-replicated fixture + /usr/bin/time -v (memory gate; cap it, see §3)
make record        # copy the record from OUT_ROOT into ./output/ (worktree runs; from the
                   # main checkout OUT_ROOT already is ./output)
```

`make test` builds the ticket-#4 synthetic fixture itself. Tests are stdlib `unittest` and
are run file-by-file in a constrained session (`-p 'test_memory.py'` etc.) per
`docs/agents/execution-environment.md`.

## 5. Real record (produced 2026-11-02)

The pinned GPU env has no parquet engine, so the derived CSV cache is built once with the
`arc` env (done; the loader fails closed with the fix command otherwise):

```bash
cd /home/chang/arc/stage2/06_prnet_transfer/06.8_baselines/execution
make tables        # arc env: 06.4_pipeline/output/input_cache from the frozen parquet tables
free -g            # abort if available < 10 GiB
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=10G -- \
  /usr/bin/time -v make run
```

Result: exit 0, **5 folds**, 951 test rows / 347 scored test conditions over the folds,
3,804 scored rows (951 × 4 baselines); all pairing invariants true; wall **5 m 37 s**
(321.6 s in-run, including hashing the 16.4 GB frozen input — 108 s system time),
**peak RSS 1,081 MB**. `make check` beforehand: 881 MiB / 21.3 s, confirming the #4
header-read loader fix.

Headline direction AUPRC@200 (per-condition macro, mean over the five folds; the primary
`no_change`/`delta_mean` comparison only — #9 does the paired inference):

| baseline | role | AUPRC@200 | up | down | RMSE | Pearson | Spearman |
|---|---|---|---|---|---|---|---|
| `no_change` | primary | 0.011011 (= 200/18,164, chance for all-tied scores) | 0.011011 | 0.011011 | 0.146374 | — | — |
| `delta_mean` | primary | **0.164489** | 0.157242 | 0.171736 | 0.141931 | 0.247906 | 0.198945 |
| `delta_mean_dose` | auxiliary | 0.160778 | 0.150774 | 0.170783 | 0.142387 | 0.240273 | 0.188591 |
| `ridge_fp` | auxiliary | 0.151178 | 0.140513 | 0.161844 | 0.150953 | 0.217464 | 0.166430 |

Per-fold `delta_mean` AUPRC@200: 0.182553 / 0.174654 / 0.161914 / 0.140450 / 0.162874.
`no_change` is exactly `k/n_genes` in every fold, so every `no_change` condition has
`pred_constant=True` and `regression_nan_reason=constant_prediction`; the fittable baselines
are non-constant with finite regression metrics. Test-row counts per fold: 192 / 188 / 191 /
190 / 190, conditions 68 / 67 / 71 / 71 / 70, disjoint from the fold's training rows, every
condition scored exactly once, and `c*`/test population identical across the four baselines.

The evaluator is reused verbatim: every `folds/fold<k>/<baseline>/eval/manifest.json` is a
`ptlib.run_manifest/v1` `kind=evaluate` record with truth
`CM per-plate matched DMSO mean (evaluation-only)` (138 control wells, 16 plates), metric
`direction_auprc` (headline k=200), and its `pred_expr_sha256` matches the prediction
directory for all 20 runs. Identity: all 13 frozen artifact sha256 (8 required re-hashed,
16,426,331,885 bytes), frozen-manifest + locked-config sha256, one gene-axis sha256 across
all folds, git `e0e2d68` (`dirty=true` because the record directory is the run's own output).

## 6. Shared-loader dependency (resolved 2026-11-02)

`ptlib.data.open_frozen_input` originally opened both h5ads through the old
`h5io.open_obs_var` (`anndata.read_h5ad(backed="r")`), which materialises every layer — the
~8.1 GB trap on the 16.2 GB Tahoe view (see `docs/agents/execution-environment.md`). Ticket
#4 replaced it with h5py/`read_elem` header reads (0.90 s / 484 MiB on that view,
regression-tested); integration `e0e2d68` carries the fix and the real `make check` above
measured 881 MiB. This step never forks the loader.
