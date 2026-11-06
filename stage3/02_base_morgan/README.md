# 02_base_morgan — structure and dose baseline

This step consumes the frozen [01_corpus](../01_corpus/README.md) release under
[ADR 0045](../../docs/adr/0045-stage3-morgan-direction-baseline.md). The current
tracer bullet implements cold-drug preparation and individual run training,
evaluation and verification. Cold-dose, warm-start and the complete 60-fit
benchmark remain subsequent work; smoke runs are explicitly marked.

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

The split seed is 20261105. Three outer folds test each parent once; validation
holds out the ceiling of 20% of each remaining parent pool (at least one parent).
All aliases and doses of a parent receive the same role across conditions.
Preparation fails when any condition–phenotype target lacks a role. Seeds 17,
42 and 2026 share these partitions. Feature metadata freezes radius 2, 2048
binary bits, chirality, effective generator options and RDKit version.

Each run saves a checkpoint, training-only log-dose statistics, learning history,
run metadata, test measurement probabilities and regenerated metrics. Evaluation
uses saved artifacts without retraining. Verification checks artifact consistency
and can bind the run to the supplied prepared release. Predictions describe
magnitude direction relative to the source's assumed same-arm vehicle=100;
class 0 is within-band, not demonstrated equivalence, and direction is not
therapeutic benefit. Seed variation is not biological replication.

The public subprocess acceptance tests use small synthetic data with aliases,
both conditions, four doses and all three classes:

```bash
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python -m unittest discover -s stage3/02_base_morgan/execution/tests -p test_workflow.py -v
```
