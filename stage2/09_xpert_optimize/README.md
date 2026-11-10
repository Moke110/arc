# Step 09 XPert gene-axis and compression ablation

[ADR 0057](../../docs/adr/0057-step09-xpert-landmark-reference.md) defines the
cold-cm-drug ablation. Round 1 screens fold 0 and seed 17. The landmark
reference follows the author’s 978-gene order. The other arms test gene-axis
size and token-axis compression against their designated controls.

## Landmark release

`execution/01_corpus/build.py` builds
`output/corpus/tahoe_tenv_bulk_l1000/` from the pinned Step 02 native count
corpora, the author’s landmark metadata and the frozen Step 03 tenv release.
The builder reproduces every train, validation and test membership in each of
the three tenv folds. Tahoe uses its eligible Step 04 profiles and matched
plate-and-cell-line DMSO controls. A SHA-256 ordering of Tahoe plate identifiers
with seed 20261110 assigns the first fifth, rounded down with a minimum of one,
to validation. All other eligible plates go to training.

Each raw well is projected to the model’s axis before CP10K and natural log1p.
Zero libraries stay zero. Controls average normalized wells within a matched
stratum, with same-plate DMSO for tenv. The identity lookup uses 256 float32
channels and seed 20261111. The HDGFL3 vector follows the CM vectors in the
same random stream and is stored separately.

The predictor includes controls, Morgan fingerprints, case metadata and
training and validation treated expression. The evaluator contains held-out
treated expression, the test case-to-native-row join, and the exact control
row joins. Predictor entry points load only their fixed, hash-verified paths.

## Landmark qualification and replay

Run from the repository root with the shared `arc` Python environment. Set
`RECORD` and `RUN` to fresh paths. The qualification policy uses a bounded
update and validation budget.

```bash
STEP=stage2/09_xpert_optimize
PY="$HOME/miniconda3/envs/arc/bin/python"
RUN="$STEP/execution/work/landmark-qualification"
RECORD="$RUN-resources"
PREDICTOR="$STEP/output/corpus/tahoe_tenv_bulk_l1000/predictor"
POLICY="$STEP/execution/landmark_record/qualification-policy.json"
CAP="$STEP/execution/script/run_capped.sh"
CLI="$STEP/execution/script/xpert.py"

"$CAP" "$RECORD/build" "$PY" "$STEP/execution/01_corpus/build.py"
"$CAP" "$RECORD/fit" "$PY" "$CLI" fit --predictor "$PREDICTOR" \
    --fold 0 --seed 17 --policy "$POLICY" --output "$RUN"
"$CAP" "$RECORD/predict" "$PY" "$CLI" predict --predictor "$PREDICTOR" \
    --checkpoint "$RUN/tenv.pt" --role validation --output "$RUN/validation"
"$CAP" "$RECORD/replay" "$PY" "$CLI" replay --predictor "$PREDICTOR" \
    --checkpoint "$RUN/tenv.pt" --prediction "$RUN/validation"
"$PY" -m unittest discover -s "$STEP/execution/tests" -p test_landmark.py -v
```

Training starts from random weights. It selects the Tahoe checkpoint by
Tahoe validation objective, initializes tenv fine-tuning from that checkpoint,
and selects the final checkpoint by tenv validation objective. The policy
fixes both update schedules, effective batch size, precision and attention
kernel. The last validated update contains optimizer and random-generator
state for a compatible `fit --resume`. Replay rebuilds per-well Δ predictions
from the selected checkpoint and compares the artifact hashes.

[Landmark evidence](execution/landmark_record/) contains the input manifests,
qualification policy and measured resource and replay records.
