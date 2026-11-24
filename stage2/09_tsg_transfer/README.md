# TranSiGen Tahoe-to-tenv transfer

See the [historical evidence index](../03_data_train/03.1_hvg-deg-5k/historical/README.md) for sealed records.

[Path migration and sealed evidence](../07_tsg_dose/log/step-path-migration.md).

[Specification #255](https://github.com/Moke110/arc/issues/255) evaluates
Tahoe pretraining and tenv fine-tuning across three tenv cold-drug-dose
folds. The study compares tenv-only training, direct Tahoe-only prediction
and Tahoe-to-tenv fine-tuning on exactly the same test profiles.

## Accepted results

See the [paired fold scores and equal-fold comparison](log/summary.md),
[comparison figure](figure/transfer.png) and
[formal acceptance receipt](log/acceptance.json).

## Protocol

All paths reuse the Step 04 32-dimensional Latent embedding and axis-local
`G_5000` expression from [ADR 0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md).
Tahoe keeps replicates of each parent, cell line and dose condition together.
Each validation parent, cell line and dose has independent training coverage.
Three Tahoe fits restore the strict minimum treated-expression validation MSE.
The lowest fold MSE selects one immutable checkpoint, with equal MSE broken
by the earlier fold ID. Tahoe selection uses Tahoe data only.

The tenv folds reuse the Step 04 cold-drug-dose split with held doses 0.2,
1 and 10 µM. The held parent × held dose intersection is the test set.
Both off-intersection held strips are excluded from fit and validation.
Both tenv-trained paths use a fresh Adam optimizer at 1e−5, the same train-only
log-dose standardizer and identical validation settings. Tahoe uses Adam at
1e−4. All fits use batch 64, training seed 42, strict best-validation
treated-expression MSE, patience 50 and at most 1000 epochs. The transfer
initialization converts the first dose affine layer to keep the original
dose embedding across physical doses before fine-tuning all model parameters.
Direct Tahoe-only predictions retain the checkpoint's original transform.

Score each test profile by gene-wise Pearson correlation between predicted
and observed Delta relative to its matched control. Give defined scores equal
weight within each fold and give all three folds equal weight in the final row.
Report undefined PCC counts alongside scores. Transfer − tenv-only is the
effect of the initialization under this fixed protocol. Transfer − Tahoe-only
is the effect of tenv adaptation from the same Tahoe checkpoint.

## Raw-first inputs and assessment

Use `execution/raw_first.py prepare` for study input freezes. It opens the
archived well-level full view via the dose step and verifies both sources
in one pair. Pass `--pair` to reuse an existing raw-first pair or `--view`
to prepare from an explicit raw-first view release. The transfer receipt pins
the pair and its view, axis and corpus digests, both representation manifests,
the shared seed values and both split tables. Tahoe validation uses the
original condition-grouped split computed from Tahoe only.

```bash
STEP=stage2/09_tsg_transfer
python "$STEP/execution/raw_first.py" prepare --candidate "$PREPARED"
python "$STEP/execution/raw_first.py" evaluate --prepared "$PREPARED" \
  --fold 0 --path transfer --prediction "$PREDICTIONS" \
  --prediction-ids "$PROFILE_IDS" --output "$ASSESSMENT"
python "$STEP/execution/raw_first.py" compare --output "$COMPARISON" "${ASSESSMENTS[@]}"
python "$STEP/execution/raw_first.py" render --data "$COMPARISON" --output "$FIGURE"
```

For an existing pair, supply `--pair "$PAIR"` to both preparation and
evaluation. Prediction arrays use treated-expression float32 values in the
order of `PROFILE_IDS`, a one-column CSV headed `profile_id`. The evaluator
uses shared metrics with the matched-control layer. Pass the nine assessment
directories for the three folds and three paths as `ASSESSMENTS`. Comparison
retains sample-equal test summaries and equal-fold means. The renderer reads
only the two exported panel CSV files through shared plotting and palette APIs.

Targeted fixture test from the repository root:

```bash
python -m unittest discover -s "$STEP/execution/tests" -p test_raw_first_transfer.py
```

## Report acceptance

Set `REPORT` to a fresh staging directory. Run each command in its own
process from the repository root.

```bash
for action in prepare verify; do
  python "$STEP/execution/publish.py" "$action" --inputs "$PAIR" \
    --native-source "$SOURCE" --candidate "$CANDIDATE" \
    --phases "$EVIDENCE/phases" --report "$REPORT"
done
python "$STEP/figure/script/render.py" --data "$REPORT/figure/data" --output /tmp/transfer-replay.png
cmp "$REPORT/figure/transfer.png" /tmp/transfer-replay.png
```

The report gate verifies nine fit/replay receipts, nine test/replay receipts,
the frozen test identities across paths, source-only selection and formal
resource evidence. The plot renderer reads only the two panel CSVs. Run
`script/check_paths.sh`, `script/check_tracked.sh` and
`script/verify_all.sh` on the data-owning checkout before release.
Run the offline manifest scan outside formal training under the Step 07 cap
and record its resource use.

```bash
stage2/07_tsg_dose/execution/capped.sh /usr/bin/time -v bash script/verify_all.sh
```

Preserve the candidate and phase evidence under the ignored recovery path
`execution/checkpoints/accepted/`. Preserve the exact five generated report
files under a dedicated report archive there. For verification from tracked
files, materialize just those five files in an empty directory with their
original layout. Run `publish.py verify --report` against that archive so the
inventory excludes scripts and execution files in the working tree.

After all gates pass, copy the staged files into their matching paths under
this step: `figure/data/transfer_scores.csv`,
`figure/data/transfer_differences.csv`, `figure/transfer.png`,
`log/summary.md` and `log/acceptance.json`.
