# Publish and replay the Step 06 diagnosis

The publication consumes terminal screening and replication decisions and
the source audit. Use the arc interpreter and run from the repository root.
Submit these commands to the Step 06 scheduler with the scoring memory cap.
The publisher streams hashes of local prediction and checkpoint binaries.
Replay also computes AP@100 from the sealed predictions and frozen truth
through the original scoring entry points, then verifies the regenerated
record, tables and figures byte for byte. It needs the same Matplotlib version and
fonts for byte-identical PNGs.

```bash
step="$PWD/stage2/06_pdiff_transfer_optimize"
cli="$step/execution/script/publication.py"
out="$step/output"
record="$step/execution/publication_record"

python "$cli" publish "$out" --record "$record"
python "$cli" replay "$out" --record "$record"
```

The step root holds the paired per-run and per-dose tables. These include
directional AP@100, checkpoint steps, valid-group counts, gains, original and
residual CM gaps, and recovery where the matched original gap is positive.
The compact JSON record retains terminal seed dispositions, local artifact paths and SHA-256
digests, and the scheduler records available before publication. The
source audit exposure table is copied to the step root. The [audit
entry point](source_audit.md) reproduces the full bracket and sensitivity
panels separately. Prediction arrays, model binaries, and scheduler
telemetry remain in the local output tree.

The attribution decision uses the primary mode only. A major contribution
requires a seed 17 gain of at least +0.010, positive gains for both added
screening seeds, a positive equal-seed mean paired effect in fold 1, and
recovery of at least half a positive original CM gap on the matched
confirmation comparison. The record includes the paired and CM-matched
seed counts for that comparison when it runs.
