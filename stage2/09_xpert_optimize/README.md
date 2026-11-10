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

## Gene-axis release

`execution/01_corpus/gene_axes.py` reproduces the Step 09 cohort from the
landmark manifest and freezes panels for all three tenv folds. It projects
eligible Tahoe training wells and the fold-local tenv training wells onto
the CM axis. Each source ranks genes by Seurat normalized dispersion of
linear CP10K, using 20 equal-width bins over log mean expression. The sum of
the two ranks gives both sources equal weight. Frozen CM order breaks ties
and orders the selected genes within each panel.

The release materializes the four model views for fold 0. Each view uses the
same rows, controls, conditioning and identity lookup. Expression is
normalized after raw-count projection to that view. The shared panel is the
intersection of the landmark and smallest G axis. Its view contains tenv
expression for scoring. The frozen shared-panel sizes are 73, 72 and 71
genes for folds 0, 1 and 2. Continuous metrics on these panels are provisional.

Use the landmark commands above with the new release and its view predictor
path. The fit CLI derives the uncompressed model from the predictor manifest.
The qualification-only prediction limit keeps GPU checkpoint replays bounded.

```bash
AXES="$STEP/output/corpus/gene_axes"
POLICY="$STEP/execution/gene_axes_record/qualification-policy.json"
"$CAP" "$RECORD/gene-axes-build" "$PY" "$STEP/execution/01_corpus/gene_axes.py"
PREDICTOR="$AXES/views/G18164/predictor"
"$CAP" "$RECORD/full-axis-fit" "$PY" "$CLI" fit --predictor "$PREDICTOR" \
    --fold 0 --seed 17 --policy "$POLICY" --output "$RUN-full"
"$CAP" "$RECORD/full-axis-predict" "$PY" "$CLI" predict --predictor "$PREDICTOR" \
    --checkpoint "$RUN-full/tenv.pt" --role validation --limit 4 --output "$RUN-full/validation"
"$CAP" "$RECORD/full-axis-replay" "$PY" "$CLI" replay --predictor "$PREDICTOR" \
    --checkpoint "$RUN-full/tenv.pt" --prediction "$RUN-full/validation"
"$PY" -m unittest discover -s "$STEP/execution/tests" -p test_gene_axes.py -v
```

[Gene-axis evidence](execution/gene_axes_record/) records the frozen views,
GPU qualification policy and checkpoint replays.

## Matrix MLP compression

All four compression arms read the `G18164` predictor. The encoder maps
18,164 gene tokens to the chosen latent width once, before the treatment
and baseline branches. The same decoder restores 18,164 tokens for each
branch before the three prediction heads. Each map has two linear layers
with GELU between them, a hidden token count equal to the latent width,
and weights shared across all 256 channels. The shared XPert parameters
initialize before the compression maps, so seed 17 gives the same initial
shared weights as `XPert_G18164`.

| Model | Latent tokens | Compression parameters | Total trainable parameters |
|---|---:|---:|---:|
| `XPert_MLP512` | 512 | 19,143,924 | 24,731,511 |
| `XPert_MLP1024` | 1,024 | 39,318,260 | 44,905,847 |
| `XPert_MLP2048` | 2,048 | 82,812,660 | 88,400,247 |
| `XPert_MLP4096` | 4,096 | 182,384,372 | 187,971,959 |

Use the full-axis predictor from the previous section. Choose one
compression arm per capped fit. Prediction and replay restore the arm
from its selected checkpoint.

```bash
ARM=XPert_MLP512
PREDICTOR="$AXES/views/G18164/predictor"
POLICY="$STEP/execution/compression_record/qualification-policy.json"
"$CAP" "$RECORD/$ARM-fit" "$PY" "$CLI" fit --model "$ARM" --predictor "$PREDICTOR" \
    --fold 0 --seed 17 --policy "$POLICY" --output "$RUN-$ARM"
"$CAP" "$RECORD/$ARM-predict" "$PY" "$CLI" predict --predictor "$PREDICTOR" \
    --checkpoint "$RUN-$ARM/tenv.pt" --role validation --limit 4 --output "$RUN-$ARM/validation"
"$CAP" "$RECORD/$ARM-replay" "$PY" "$CLI" replay --predictor "$PREDICTOR" \
    --checkpoint "$RUN-$ARM/tenv.pt" --prediction "$RUN-$ARM/validation"
"$CAP" "$RECORD/compression-fixture" "$PY" -m unittest discover \
    -s "$STEP/execution/tests" -p test_compression.py -v
```

[Compression evidence](execution/compression_record/) records all four
real GPU qualification fits and checkpoint replays.
