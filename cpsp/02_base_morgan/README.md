# 02_base_morgan — structure and dose baseline

This CarPS-Phenotype step establishes a structure-and-dose comparator for
cardiomyocyte viability direction. It fits separate models for the baseline and
hypoxia–reoxygenation (H-R) arms of the frozen Heliste release. The models
estimate the probabilities of decrease, within-band response and increase
under the corpus's same-arm vehicle assumption. These classes describe assay
response magnitude. [Heliste corpus](../01_corpus/README.md) ·
[Morgan protocol](../../docs/adr/0045-stage3-morgan-direction-baseline.md)

## Prepare the release and fit the baseline

Preparation verifies the corpus release and retains each admitted source
measurement, including aliases. Each model takes a chirality-aware, radius-2,
2,048-bit Morgan fingerprint of the desalted parent and the log dose scaled
using only training measurements. A 128-unit ReLU hidden layer feeds three
class logits in the order `[-1, 0, +1]`. Validation macro-F1 selects the
checkpoint. [Corpus manifest](../01_corpus/output/manifest.json) ·
[Training implementation](execution/training.py)

The split seed is `20261105`. Seeds `17`, `42` and `2026` reuse the same
partitions, with identical group roles in both arms.

| Setting | Folds | Test partition |
| --- | ---: | --- |
| Cold drug | 3 | Holds out entire parents, including all doses and aliases. |
| Cold dose | 4 | Holds out one dose globally and retains every test parent in training. The 0.04 and 5 µM folds test extrapolation. The 0.2 and 1 µM folds test interpolation. |
| Warm start | 3 | Holds out parent–dose combinations while retaining each parent and dose in training. |

The [split manifest](execution/split_manifest.json) records all ten partitions
over 5,488 measurements. The two arms and three training seeds give 60 formal
fits.

## Formal results and release evidence

All 60 fits completed, with zero failed or pending. Cold drug had the lowest
mean macro-F1 in both arms. The +1 class had zero recall in all cold-drug
fits. Values below are unweighted fold-by-seed means ± sample SD. The installed
report at `output/formal/report.md` also gives each held-dose result, log loss,
classwise precision and recall, and pooled cold-drug results by seed.
[Release record](execution/release_record.json)

| Setting | Baseline macro-F1 | H-R macro-F1 |
| --- | ---: | ---: |
| Cold drug | 0.4065 ± 0.0229 | 0.4115 ± 0.0394 |
| Cold dose | 0.5635 ± 0.0219 | 0.5623 ± 0.0449 |
| Warm start | 0.5504 ± 0.0448 | 0.5301 ± 0.0149 |

The independent audit reproduced all run metrics and aggregates from 49,392
saved test predictions. The original formal execution passed both `evaluate`
and `verify`. The CarPS-Phenotype migration preserved the split manifest,
formal report and all 60 runs' artifacts. On the migrated release,
`verify --require-complete` replayed all 60 saved fits and regenerated the
summary without refitting.
[Original execution](execution/validation/formal-execution.md) ·
[Independent audit](execution/validation/formal-data-audit.md) ·
[Replay implementation](execution/benchmark.py)

## Downstream comparison

A future CarPS-Phenotype experiment can compare the Morgan baseline with a
predictor that takes a CarPS-Transcriptome predicted response Δ and dose.
The Morgan release provides the reference partitions, labels, metric
conventions and per-measurement probabilities for paired evaluation on rows
covered by both predictors. Integration and Heliste coverage of transcriptome
predictions remain pending. [Comparison design](../../docs/adr/0045-stage3-morgan-direction-baseline.md) ·
[CarPS-Transcriptome](../../cpst/README.md)

## Run and verify

Use the step instructions for preparation, fitting, formal replay, synthetic
acceptance tests and installed-release verification. [Step instructions](AGENTS.md#run-and-verify)
