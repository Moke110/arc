# Formal Heliste Morgan baseline execution (#63)

The frozen matrix contains 60 completed formal fits and zero failed fits: two
condition–phenotype targets × three training seeds × (three cold-drug folds,
four cold-dose folds, three warm-start folds). All fits use the frozen schedule;
there was no hyperparameter search, smoke-epoch override or biological replication.

## Frozen release and interpretation

`../release_record.json` records every run identity, all five persisted artifact
SHA256 values, selected epoch, training-only dose statistics, execution versions,
config, feature generator options, source fingerprints and metric conventions.
`../split_manifest.json` persists the actual ten partitions and their parent or
parent–dose assignments, plus the 5,488 ordered measurement IDs, condition, dose,
parent, label and source label. Each split's roles array indexes that measurement
order. Seeds 17, 42 and 2026 reuse split seed 20261105.

The authenticated frozen Heliste release admits 2,744 measurements per condition,
686 resolved source labels and 663 canonical parent structures; 12 source values
per condition remain pending and are not silently imputed. Class order is
[-1, 0, 1]: magnitude below 80, inclusive [80,120], above 120 under the explicit
same-arm vehicle=100 assumption. Class 0 does not establish equivalence; direction
does not establish therapeutic benefit. The source values lack the evidence to
claim statistical significance. Seed variation is not biological replication.

Features are 2,048 binary radius-2 Morgan bits with chirality and log10 dose,
with training-only measurement-weighted population standardization (ddof=0).
The model has one 128-unit ReLU hidden layer and three output logits, trained
with AdamW and unweighted cross entropy on logits; batch size 64, learning rate
0.001, weight decay
0.0001, maximum 300 epochs and patience 30. Saved validation selection and test
predictions are preserved. Metric labels, absent-class handling and float64
log-loss epsilon are frozen in the release record. Aggregates use sample standard
deviation (ddof=1); a single observation has null standard deviation.

## Commands and resource accounting

Executed from the repository root with the `arc` Python executable:

```bash
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
export CUDA_VISIBLE_DEVICES=''
python stage3/02_base_morgan/execution/workflow.py prepare --input stage3/01_corpus/output --output stage3/02_base_morgan/output/prepared --settings cold_drug cold_dose warm_start
python stage3/02_base_morgan/execution/workflow.py benchmark --prepared stage3/02_base_morgan/output/prepared/prepared.json --output stage3/02_base_morgan/output/formal
```

The production runner serializes CPU fits, with two library threads, CUDA hidden,
a hard 10 GiB RLIMIT_AS child cap and five-second monitoring. Host-visible
telemetry enforces available RAM ≥12 GiB and excludes concurrent external trainers.
The run paused for external Stage 2 probes and resumed the same pending-only
ledger during quiet windows. Maintenance and source-only validation changes did
not trigger additional fits. No external process was terminated. Completed
artifacts were verified and reused. Initial code hashes and explicit ae2426a
(completion/thread handling) and a33ffe1 (resource classification) transitions
are preserved in the release record; training.py, frozen input, split, features
and model schedule stayed unchanged. The final pending-only driver exited 0.

Ignored `output/formal/resources/*/{resources.jsonl,time.txt}` records per-fit
telemetry and GNU time evidence. Maximum single-fit GNU time RSS was
889,835,520 bytes (848.61 MiB); maximum sampled child-tree RSS was 777,007,104
bytes. The orchestrator retains its own Python/Torch memory: summing the
workflow driver and recursive PID descendants in each of 240 raw host process
snapshots gives a maximum 1,552,723,968 bytes (1.446 GiB). This sampled maximum
is not a kernel peak and is distinct from the per-fit peak. Minimum sampled
available host RAM was 43,445,182,464 bytes. Maximum sampled host GPU occupancy
was 1,823 MiB; CPU fits have CUDA hidden and host occupancy includes desktop or
external activity. Paused external probes are documented in the execution logs,
not represented as model GPU use.

## Public replay

Evaluation and verification use saved artifacts without training, with the same
thread/CUDA settings and `/usr/bin/time -v`. Each shell applies `ulimit -v
10485760` before executing the command:

```bash
python stage3/02_base_morgan/execution/workflow.py evaluate --benchmark stage3/02_base_morgan/output/formal --require-complete
python stage3/02_base_morgan/execution/workflow.py verify --benchmark stage3/02_base_morgan/output/formal --prepared stage3/02_base_morgan/output/prepared/prepared.json --require-complete
```

Both public replay commands exited 0. Evaluation took 3:26.68
with maximum RSS 856,592 KiB; verification took
3:28.07 with maximum RSS 856,260 KiB.
All 300 persisted files from the 60 completed fits remain byte-identical to the
pre-replay snapshot, including all files in the original 19-fit and 28-fit
snapshots. Summary and report SHA256 values remain unchanged after verification.
The release record binds replay logs, GNU time files and hash snapshots. Raw
command output lives in ignored `output/public-{evaluate,verify}.log`; GNU time
logs live in `output/public-{evaluate,verify}-time.txt`. Full metrics and
setting/seed aggregates live in ignored `output/formal/{summary.json,report.md}`.
Independent data audit at commit fbbfaf1 checked all 49,392 predictions, per-run
metrics and aggregates successfully. No replay retrained a model.
