# Independent formal data and numerical audit — issue #63

Audited the actual frozen corpus, prepared release and all 60 completed formal
runs on `dev`, 2026-11-06. The auditor used a standalone temporary Python script
with standard-library CSV/JSON/hash routines, NumPy and RDKit. It imported no
feature preparation, split, training, reporting or metric helpers, and no Torch.
No model was trained or altered by this audit.

## Data and partitions

- All eight files named by `stage3/01_corpus/output/manifest.json` matched their
  streamed SHA256 checksums. Every field of all 5,488 CSV measurements matched
  the corresponding prepared row by condition, source asset and source cell;
  no measurement, alias, readout, label or provenance was lost or averaged.
- Each condition retained 2,744 measurements, 686 source labels and 663 parents.
  Prepared measurement IDs were unique. RDKit independently regenerated all 663
  parent InChIKeys from their canonical SMILES. Reapplying the frozen `<80`,
  `[80,120]`, `>120` rule reproduced every direction label.
- The canonical JSON SHA256 of the complete prepared object matched the ledger's
  prepared fingerprint. The split seed was 20261105; training seeds were
  17, 42 and 2026. The prepared and ledger identities exactly equaled the
  60-member Cartesian protocol: two condition–phenotype targets times three
  seeds times (three cold-drug, four cold-dose, three warm-start folds).
- Every manifest assigned every measurement exactly one role, kept aliases
  together, shared group roles across conditions, and reproduced all persisted
  role measurement/parent/combination/class counts. Each target had all roles.
  Each measurement was tested once across the outer folds of each setting.
- Cold-drug parent roles were disjoint. Cold-dose held levels were globally
  absent from training/validation; held levels and interpolation/extrapolation
  tags matched 0.04, 0.2, 1 and 5 µM. Every validation/test parent remained in
  training for cold-dose and warm-start. Warm-start also retained every held
  dose marginal in training, with disjoint parent–dose combinations.

Group counts below are train / validation / test; membership is shared across
conditions. Every validation count equals the ceiling of 20% of non-test groups.

| Setting | Folds | Group | Counts per fold |
|---|---:|---|---|
| Cold drug | 3 | Parent | 353 / 89 / 221 |
| Cold dose | 4 | Parent–dose | 1,591 / 398 / 663 |
| Warm start | 3 | Parent–dose | 1,414 / 354 / 884 |

## Saved run evidence

All 60 run manifests reported completed formal execution with `smoke_epochs`
null. Their identity, input/split fingerprints and feature metadata matched the
prepared release. All 240 artifact checksums (checkpoint, history, predictions,
metrics per run) matched streamed hashes.

Across 49,392 saved test predictions, measurement IDs exactly matched the
independently selected test rows with no duplicates; every truth, source field
and run identity matched. All probability vectors had three finite entries in
[0,1], summed to one within absolute tolerance 1e-6, and produced the saved hard
class by first-index argmax in [-1,0,+1] order.

For every run the auditor independently rebuilt the confusion matrix, per-class
precision/recall/F1 and fixed-three-class macro-F1, using zero for zero
denominators. It recomputed log loss in float64 after clipping to machine epsilon
and 1−epsilon and normalizing each probability row. All saved metrics matched
within 1e-12 absolute and relative tolerance.

Training dose means and population standard deviations were independently
recomputed from only that target's training measurement rows; counts, ddof=0 and
saved values matched. Every configuration exactly matched 2049→128→3, AdamW
learning rate 0.001, weight decay 0.0001, batch size 64, maximum 300 epochs and
patience 30. Histories contained sequential epochs, finite losses/scores and
31–87 epochs. Every selected epoch was the earliest validation macro-F1 maximum;
every history stopped exactly after 30 unimproved epochs (or the 300-epoch cap),
with no earlier patience violation. This verifies persisted selection evidence;
architecture/checkpoint inference is covered by the separate workflow/code review.

## Aggregate audit

The final summary reported 60 completed formal runs and complete protocol scope.
The independent audit reproduced every run metric and all 14 aggregate records:
unweighted overall fold–seed means, separate held-dose means, and their sample
standard deviations (`ddof=1`) for confusion matrices, precision, recall, F1,
macro-F1 and log loss. Each contributing run count matched. It independently
concatenated the three cold-drug outer tests for each target and seed, checking
six unique 2,744-measurement out-of-fold sets and both across-seed aggregates.
All numerical comparisons passed at 1e-12 absolute/relative tolerance. Pooled
measurement metrics were checked independently of the unweighted fold means.

Audited final `summary.json` SHA256:
`1ea09cba4d6b186d9ebf4df7bc304c84b39caf2c0cc8007a0625c072d14fdd4c`.
The stale intermediate summary was explicitly rejected as completion evidence.

## Resources and command

`free -g` reported 40 GiB available before loading. The preferred systemd scope
could not connect to the user bus in this sandbox, so the documented fallback
used a 10 GiB virtual-memory cap and two threads for each numerical library:

```bash
bash -c 'ulimit -v 10485760; exec /usr/bin/time -v env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 /home/chang/miniconda3/envs/arc/bin/python /tmp/formal_data_audit.py'
```

The final all-artifact and aggregate audit exited 0 in 0.76 s, peaking at
138,520 KiB RSS with zero swaps. The earlier per-run pass used 138,560 KiB RSS
and explicitly waited for summary regeneration. The temporary auditor and machine-readable
results reside under `/tmp`; no feature code, source data or formal artifacts
were modified. No data, partition, artifact or numerical discrepancy was found.
