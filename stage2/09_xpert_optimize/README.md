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

## Round 1 evaluation and replay

The [frozen plan](execution/round1_record/plan.json) pins the cohort, fit policy,
and all nine model entries. The two deterministic references use the full CM
axis and the shared panel. `delta_mean` first averages training wells within
each tenv drug–dose group, then averages the groups equally. Both references
use zero training minutes.

The evaluator ranks genes by absolute predicted group Δ. Equal scores enter
together. `AUPRC@full` ranks covered CM genes against all measured significant
CM genes. `AUPRC@G` ranks the model axis against significant genes on that axis,
including HDGFL3 for the landmark reference. A measured, testable group with
zero significant genes on the relevant axis has undefined AP. Each PR curve
averages precision over groups that reach each recall and shows the count of
contributing groups. The reported macro AP averages per-group AP.

Continuous metrics compare group Δ on the shared panel. Per predicted well,
the evaluator adds predicted native-axis Δ to its matched native-axis control,
converts to nonnegative expression, and projects and normalizes on the
shared panel. It subtracts the control independently normalized on that
panel. Observed shared-panel Δ comes from independently projected raw counts.
Native-axis continuous scores are reported separately in the tables.

`build_significance.py` builds the two-sided DESeq2 label releases with
same-plate DMSO and within-contrast BH correction over CM plus HDGFL3. The
validation seal binds the replayed validation predictions, labels and scores.
It records each strict `AUPRC@full` improvement against the designated
control. Test labels and scoring require the seal committed on `dev`.

```bash
ROUND="$STEP/output/round1"
WORK="$STEP/execution/tmp/round1"
PROMOTION="$STEP/execution/round1_record/promotion.json"
RSCRIPT="$HOME/miniconda3/envs/biomni_e1/bin/Rscript"
QUEUE="$STEP/execution/script/round1_fit.py"
SCORE="$STEP/execution/script/evaluate.py"
LABELS="$STEP/execution/script/build_significance.py"

"$PY" "$QUEUE" --role validation --work-root "$WORK"
"$CAP" "$RECORD/validation-labels" "$PY" "$LABELS" --role validation \
    --rscript "$RSCRIPT" --output "$ROUND/labels/validation"
"$CAP" "$RECORD/validation-score" "$PY" "$SCORE" seal --work-root "$WORK"
# Commit the validation promotion on dev before running the test phase.
"$PY" "$QUEUE" --role test --work-root "$WORK"
"$CAP" "$RECORD/test-labels" "$PY" "$LABELS" --role test \
    --rscript "$RSCRIPT" --workers 2 --promotion "$PROMOTION" --output "$ROUND/labels/test"
"$CAP" "$RECORD/test-score" "$PY" "$SCORE" test --work-root "$WORK"
"$CAP" "$RECORD/evaluation-replay" "$PY" "$SCORE" replay --work-root "$WORK"
```

Both score releases include aggregate metrics, per-group scores and PR
data. The test release also publishes `xpert_opti_r1_metrics.png` and
`xpert_opti_r1_PR.png` at the Round 1 root. Independent replay recalculates
the validation seal, test tables and both figures from frozen inputs and
checks every output by SHA-256.

## Round 1 results

All nine XPert arms completed the formal fold 0, seed 17 policy and replayed
their selected checkpoints on tenv validation and test. The validation seal in
`88b621b` advanced all four G arms and all four MLP arms for confirmation.

On the exploratory test split, `XPert_MLP2048` led the trained arms on
AUPRC@full at 0.02823, compared with 0.02480 for `XPert_G18164` and 0.02984 for
`delta_mean`. AUPRC@full was defined for 359 of the 462 drug–dose groups.
The shared continuous panel contains 73 genes, so its metrics are provisional.

- [Metrics figure](output/round1/xpert_opti_r1_metrics.png)
- [Precision–recall figure](output/round1/xpert_opti_r1_PR.png)
- [Test scores and PR data](output/round1/test_scores/)
- [Acceptance, checkpoint replay and resources](execution/round1_record/acceptance.json)

The work directory holds fits, queue status and per-attempt resource logs.
Sealing requires a terminal validation outcome for every model. Test scoring
requires a terminal test outcome for every model whose validation replay
succeeded. A terminal failure remains visible in the score release.
The tables distinguish checkpointed training phase minutes from total known
fit-attempt wall minutes and list how many attempts have timing evidence.

## Round 1 learning curves and endpoint convergence check

Each trained arm has its own `output/round1/xpert_opti_r1_learning_XPert_*.png`,
with separate Tahoe pretraining and tenv fine-tuning panels. Blue is the
training objective averaged over the preceding validation interval; green is
the validation objective at the indicated update. Both use the frozen loss
weights `(0.2, 0.003, 0.2, 1.0)`. Open circles mark the best recorded validation
point. The deterministic references do not train and have no learning curves.

Reproduce from existing fit JSON logs, without loading data or checkpoints:

```bash
"$HOME/miniconda3/envs/arc/bin/python" stage2/09_xpert_optimize/execution/script/round1_learning_curves.py
```

[Curve points](output/round1/learning_curve_points.csv) and
[endpoint audit](output/round1/learning_curve_convergence.csv) retain the
numeric evidence. [Source hashes](execution/round1_record/learning-curves.json)
bind the original fit logs. This audit is separate from held-out score replay.
The endpoint criterion requires final validation loss to be no lower than the
previous observation and strictly above its recorded minimum.

**None of the nine arms passes this criterion on tenv:** all improve between
updates 300 and 600 and reach their recorded minimum at the final update.
On Tahoe, the four MLP arms pass the endpoint criterion (best at update 256,
higher loss at 512); the base and four G arms continue improving through 512.
Only two observations exist per phase, so these are sparse learning curves
and an endpoint check, not evidence of asymptotic convergence. More updates
and more frequent validation are needed to establish a plateau on tenv.
