# CPSP 02_base_morgan agent instructions

Read the [step README](README.md) for the formal results and release evidence. Follow the frozen Morgan protocol in [ADR 0045](../../docs/adr/0045-stage3-morgan-direction-baseline.md). Before loading corpora or fitting, activate `arc` and apply the [workstation limits](../../docs/agents/execution-environment.md).

## Run and verify

- **Prepare.** The [public workflow](execution/workflow.py) verifies the materialized `cpsp/01_corpus/output` release before reading its rows. Check the [corpus manifest](../01_corpus/output/manifest.json) for the input inventory. Use a single CSV for synthetic acceptance. Preparation defaults to cold drug. Request all three settings for the full benchmark and reuse a preparation directory only when its contents match.
- **Fit.** Use `train` with the prepared JSON for one condition–phenotype, setting, fold and seed. Use `benchmark` for the sequential, resource-monitored matrix. After diagnosing an interrupted or failed fit, rerun the same benchmark with `--retry-failed` to archive the attempt and retain verified completed runs. Reserve `--smoke-epochs` for integration checks.
- **Replay.** Use `evaluate` to regenerate reports from saved artifacts. Use `verify --benchmark` with `--prepared` and `--require-complete` to check all 60 formal runs against the prepared release. Both benchmark replay modes rewrite the ledger, summary and report. To replay a sealed release, copy its `prepared` and `formal` directories together into scratch and preserve the checkpoint symlinks.

From the repository root, prepare, fit and verify a new formal benchmark in scratch.

```bash
mkdir -p tmp
MORGAN_WORK_DIR=$(mktemp -d "$PWD/tmp/morgan.XXXXXX")
python cpsp/02_base_morgan/execution/workflow.py prepare \
  --input cpsp/01_corpus/output --output "$MORGAN_WORK_DIR/prepared" \
  --settings cold_drug cold_dose warm_start
python cpsp/02_base_morgan/execution/workflow.py benchmark \
  --prepared "$MORGAN_WORK_DIR/prepared/prepared.json" --output "$MORGAN_WORK_DIR/formal"
python cpsp/02_base_morgan/execution/workflow.py verify \
  --benchmark "$MORGAN_WORK_DIR/formal" --prepared "$MORGAN_WORK_DIR/prepared/prepared.json" \
  --require-complete
```

For CLI changes, run the synthetic [single-run acceptance tests](execution/tests/test_workflow.py) and [matrix acceptance tests](execution/tests/test_benchmark.py).

```bash
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
  python -m unittest discover -s cpsp/02_base_morgan/execution/tests -p test_workflow.py -v
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
  python -m unittest discover -s cpsp/02_base_morgan/execution/tests -p test_benchmark.py -v
```

For an installed release, run the [migration verifier](execution/migrate_release.py) to compare the materialized artifacts with the [release record](execution/release_record.json) and pinned Git history. The checks include the unchanged [split manifest](execution/split_manifest.json) and original run hashes. For public path changes, run the [CPSP path audit](../../script/check_cpsp_paths.sh).

```bash
python cpsp/02_base_morgan/execution/migrate_release.py verify
bash script/check_cpsp_paths.sh
```

## Protocol invariants

- Preserve each measurement and source alias. Group partition roles across conditions by desalted parent for cold drug and by parent–dose for cold dose and warm start. Hold the test dose out globally for cold dose while retaining parent coverage in training. Keep both the parent and dose represented in training for warm start. Every condition–phenotype target needs all three partition roles. Reuse the same splits across seeds.
- Keep one model per condition–phenotype target, taking the chirality-aware 2,048-bit Morgan fingerprint and `log10(dose_um)` standardized with training rows only. Preserve the saved class order `[-1, 0, +1]` and use the [corpus label contract](../../docs/adr/0044-stage3-source-first-corpora.md) and [owned metric implementation](execution/metrics.py) for model comparisons. Use the frozen partition identities for a future transcriptomic predictor and pair its evaluation with Morgan on jointly covered measurements.
