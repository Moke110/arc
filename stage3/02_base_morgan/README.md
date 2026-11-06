# 02_base_morgan — structure and dose baseline

This step consumes the frozen [01_corpus](../01_corpus/README.md) release under
[ADR 0045](../../docs/adr/0045-stage3-morgan-direction-baseline.md). The public workflow implements cold-drug, cold-dose and warm-start preparation
and individual or complete matrix training, evaluation and verification.
Smoke runs are explicitly marked.

```bash
python stage3/02_base_morgan/execution/workflow.py prepare --input stage3/01_corpus/output --output /tmp/morgan-prepared
python stage3/02_base_morgan/execution/workflow.py train --prepared /tmp/morgan-prepared/prepared.json --output /tmp/morgan-run --condition baseline --phenotype cell_viability --setting cold_drug --fold 0 --seed 17 --smoke-epochs 2
python stage3/02_base_morgan/execution/workflow.py evaluate --run /tmp/morgan-run
python stage3/02_base_morgan/execution/workflow.py verify --run /tmp/morgan-run --prepared /tmp/morgan-prepared/prepared.json
```

Use the `arc` environment. Before training, read the repository
[execution limits](../../docs/agents/execution-environment.md), check available
memory and cap the process and library threads. Long training runs execute
serially; the example's two-epoch smoke run is only an integration check.
Omit `--smoke-epochs` for the frozen training schedule.

Directory inputs must pass the public corpus integrity verifier before reading.
A single CSV with the corpus fields is accepted for synthetic fixtures; this does
not authenticate it as the frozen Heliste release. Preparation preserves each
measurement and alias, canonical parent structure, InChIKey identity, provenance,
label and dose, sorted by stable measurement identity. Canonical JSON SHA256
fingerprints freeze inputs and partitions. Reusing an output requires identical
preparation; choose fresh locations for incompatible artifacts.

Preparation defaults to `--settings cold_drug`; select any combination of
`cold_drug cold_dose warm_start` explicitly to freeze the desired run scope.
Settings are stored in canonical order, independently of CLI argument order.
Completing a selected scope sets `requested_scope_complete`; all-formal outcomes
within that scope also set `formal_scope_complete`. Full `complete` and
`formal_benchmark_complete` require the canonical 60 formal fits across both
current targets, all three settings and seeds. `--require-complete` enforces
that full matrix, even when the prepared scope contains only cold-drug runs.

The split seed is 20261105. Three cold-drug outer folds test each parent once; validation
holds out the ceiling of 20% of each remaining parent pool (at least one parent).
All aliases and doses of a parent receive the same role across conditions.
Cold dose creates four global held-dose folds in order 0.04, 0.2, 1 and 5 µM.
The two boundary doses are extrapolation; the interior doses are interpolation.
Validation reserves the ceiling of 20% of non-test parent–dose combinations,
keeping every validation and test parent represented in training within each
condition–phenotype target. All aliases share the combination role. The manifest
persists held dose, generalization type, combination assignments and role counts.
Missing dose levels or impossible parent coverage fail with a diagnostic.

Preparation fails when any condition–phenotype target lacks a role. Seeds 17,
42 and 2026 share these partitions. Feature metadata freezes radius 2, 2048
binary bits, chirality, effective generator options and RDKit version.

Each run saves a checkpoint, training-only log-dose statistics, learning history,
run metadata, test measurement probabilities and regenerated metrics. Dose scaling
uses measurement-weighted population statistics (`ddof=0`); a constant training
log-dose uses scale 1. The saved statistics apply to every inference row.
Evaluation uses saved artifacts without retraining. Verification checks artifact consistency
and can bind the run to the supplied prepared release. Predictions describe
magnitude direction relative to the source's assumed same-arm vehicle=100;
class 0 is within-band, not demonstrated equivalence, and direction is not
therapeutic benefit. Seed variation is not biological replication.

The public subprocess acceptance tests use small synthetic data with aliases,
both conditions, four doses and all three classes:

```bash
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python -m unittest discover -s stage3/02_base_morgan/execution/tests -p test_workflow.py -v
```

Run the complete frozen matrix sequentially (60 fits for the current two targets):

```bash
python stage3/02_base_morgan/execution/workflow.py prepare --input stage3/01_corpus/output --output /tmp/morgan-all --settings cold_drug cold_dose warm_start
python stage3/02_base_morgan/execution/workflow.py benchmark --prepared /tmp/morgan-all/prepared.json --output /tmp/morgan-benchmark
python stage3/02_base_morgan/execution/workflow.py evaluate --benchmark /tmp/morgan-benchmark --require-complete
python stage3/02_base_morgan/execution/workflow.py verify --benchmark /tmp/morgan-benchmark --prepared /tmp/morgan-all/prepared.json --require-complete
```

`benchmark.json` records every expected identity and its pending, running,
completed or failed status. Compatible completed runs are verified and reused;
incompatible scope, inputs, schedules and run artifacts are rejected. Resource
preflight blocking leaves unstarted runs pending and stops dispatch. After
investigating a failed or interrupted run, use `benchmark --retry-failed` with
the same arguments to archive its partial artifacts and resource logs before
starting a fresh attempt. Completed runs remain verified and reused; failed
attempts remain under `failed_attempts/`. Evaluation regenerates `summary.json`
and `report.md` from persisted predictions and checks checkpoint integrity without
refitting. Missing, duplicate or unexpected run identities fail verification.

The runner launches one CPU child at a time with a 10 GiB `RLIMIT_AS` hard cap,
two library threads and CUDA hidden. Every child records `/usr/bin/time -v`
maximum RSS and monitored resource snapshots. Launches require at least 12 GiB
available RAM, observable host process and GPU telemetry, and no external long
trainer. Unknown telemetry fails closed. Do not run the formal benchmark while
another training job is active. A monitoring violation stops only the runner's
owned child process group and persists failure accounting.

For synthetic acceptance only, add `--smoke-epochs 1` to `benchmark`; its artifacts
are explicitly smoke and cannot complete a formal benchmark. The full synthetic
60-fit test substitutes the child launcher in a test harness, runs real tiny
models, and labels resource measurements as simulated. Production CLI commands
have no telemetry bypass. Run the narrow acceptance test with:

```bash
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python -m unittest discover -s stage3/02_base_morgan/execution/tests -p test_benchmark.py -v
```
