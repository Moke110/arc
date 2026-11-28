# CarPS-Transcriptome workflow

Read the owning step's guide and report before changing it. Apply the
[accepted-product contract](../docs/adr/0062-stage2-accepted-products-and-reproducible-figures.md)
and the [active layout decision](../docs/adr/0105-stage-cleanup.md). Before
loading data or training, apply the [WSL execution limits](../docs/agents/execution-environment.md).
Use the numbered guides to materialize inputs and verify releases.

## Active steps

| Guide | Owns |
| --- | --- |
| [01 Raw data](01_data_raw/AGENTS.md) | Source-native integer counts, gene coverage and control matching. |
| [02 Prediction axes](02_gene_axis/AGENTS.md) | Ordered HVG-DEG and Denoise-DEG axes and their selection evidence. |
| [03 Training data](03_data_train/AGENTS.md) | Eight tenv plate-pooled releases, frozen folds and the Tahoe well-level HVG-DEG 5k release. |
| [04 TranSiGen axis comparison](04_tsg_axis/AGENTS.md) | The six-axis tenv study, run through [matrix.py](04_tsg_axis/execution/matrix.py) and [evaluate.py](04_tsg_axis/execution/evaluate.py). |
| [05 TranSiGen loss](05_tsg_loss/AGENTS.md) | Δ loss and reconstruction objective studies on Denoise-DEG 5k. |
| [06 TranSiGen modules](06_tsg_module/AGENTS.md) | Expression widths, depths, prediction paths and measured encoder targets. |
| [07 CarPS-Transcriptome loss](07_cpst_loss/AGENTS.md) | The additive loss comparison on the shared signed Δ model. |
| [08 Chemical representations](08_cpst_chem-emb/AGENTS.md) | The comparison of seven frozen representations on the fixed alpha. |
| [09 CarPS-Transcriptome axes](09_cpst_axis/AGENTS.md) | The eight-axis UniMol comparison with the study's standardized log dose. |
| [10 CarPS-Transcriptome dose](10_cpst_dose/AGENTS.md) | The cold-dose-drug comparison and selected `no-dose` configuration. |

For shared identities, frozen representations, seeds, vector metrics, derivation,
publication or plotting, use the [utility guide](util/AGENTS.md) and
[API inventory](util/README.md). Steps 04–10 call the
[shared training and evaluation API](util/README.md#transigen-training-and-evaluation)
through their numbered entry points. Keep population and fold validation,
calibration, aggregation and report tables with the study. Bind each run to
the study's pinned data, ordered axis, split, fold, seed, chemical representation,
dose transform and conditioning, architecture, objective and source. Use the
[selected model report](README.md#selected-carps-transcriptome-configuration)
for the successor configuration and its decisions.

## Step contract

- Place the minimum reusable accepted products under `output/` at stable paths
  with one owner per product. This includes corpus and gene-axis releases,
  evaluator truth, reusable checkpoints and public reader or metric code.
  Build and audit code imports the public implementation.
- Keep build and audit code, configurations, tests and disposable intermediates
  in `execution/`. Keep audit evidence in `log/`, including inputs and hashes,
  seeds, environment, resource measurements, migration mappings and compact
  failure records. Maintain a step guide for commands and a README for
  accepted results, limitations and downstream products.
- Publish one current set of PNGs under `figure/` with complete plotted
  values in `figure/data/` and CSV-only renderers in `figure/script/`.
  Record the source of each derived plot table. Render solely from those CSVs
  and replace figures and render receipts at their canonical paths.
- Draw report marks through the [shared plotting API](util/README.md#plotting).
  Plot learning curves in one panel per model configuration, with train and
  validation traces for every fold and seed. Default to bars for other
  reports. Show each defined fold and seed as a dot, the group mean as the
  bar height and sample SD when at least two values are defined. Use hue
  for the first grouping variable and depth for the optional second.
  Keep dots, complete SD intervals and relevant references in view. Use
  common limits for directly comparable panels. Include zero in every
  comparison metrics panel and label other bar charts that omit zero.
  [Plotting decision](../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md) ·
  [Comparison standard](../docs/adr/0109-cpst-seven-metric-comparison-figures.md) ·
  [Comparison CSV schema](util/README.md#comparison-figures)

## Acceptance and replacement

Build a candidate outside the stable output, then verify its structure and
data constraints, scientific acceptance, downstream readers and CSV-only
figure reproduction before replacement. Record stochastic metric tolerances
and byte hashes separately. Keep the accepted release in place if a candidate
fails. Preserve superseded tracked evidence in Git and resolve dependencies
before discarding superseded payloads. Keep external `source/` read-only and
put declared additions in `local/` under the
[external provenance contract](../docs/adr/0030-external-data-provenance-and-step-layout.md).

## Verify the current figures

Activate `arc` at the repository root. Match Python, Matplotlib, NumPy and
any other package versions in the render receipts' `environment` fields.
Run the shared gate after figure changes. It checks figure inventory,
layout and receipt hashes, and requires two independent CSV-only renders
to match all 45 published PNGs. The comparison tests check panel order,
replicate marks and sample SD endpoints against the published tables.
[Figure verifier](util/figure/script/verify_report_figures.py) ·
[Comparison tests](util/execution/tests/test_comparison.py) ·
[Published figure evidence](README.md#verify-the-current-figures)

```bash
(
  set -e
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/cpst-figure-check.XXXXXX")
  python cpst/util/figure/script/verify_report_figures.py --work "$work/replay"
  TMPDIR="$work" MPLCONFIGDIR="$work/replay/matplotlib" \
    PYTHONPATH=cpst/util/output/api \
    python -m unittest discover -s cpst/util/execution/tests -p test_comparison.py
)
```

## Verify the active inputs

Use the materialized Steps 01–03 releases, including all eight tenv views,
their frozen folds and the Tahoe release. Work from the repository root
in `arc`, with `cpst/history` absent. Define `capped_step04` using the
[Step 04 runtime guide](04_tsg_axis/AGENTS.md#prepare-the-capped-runtime).
Require `result: PASS` and `history_available: false` in the resulting
JSON. The verifier exercises matched native tenv counts, ordered
projections and both tenv and Tahoe training readers.
[Active-input verifier](04_tsg_axis/execution/verify_active_inputs.py) ·
[Archive-isolated receipt](04_tsg_axis/log/323-archive-isolated-smoke.json)

```bash
capped_step04 python cpst/04_tsg_axis/execution/verify_active_inputs.py
```

## Replay sealed studies

Keep sealed inputs read-only. For replay, copy the fit to scratch and use its
original registered source revision and protocol. Use the
[pinned Step 06.4 replay](06_tsg_module/06.4_round4/AGENTS.md#replay-a-sealed-step-064-fit)
for registrations that bind the original `stage2/` layout. Use the
[TranSiGen verification inventory](util/training/verification.md) for
Steps 04–06 and the [CarPS-Transcriptome verification inventory](util/training/cpst-verification.md)
for Steps 07–10. Retrieve tracked historical code, figures and receipts from
their original Git revisions.
[ADR 0087 cutover](https://github.com/Moke110/carps/blob/99078e6f/stage2/README.md#deploy-the-adr-0087-cutover) ·
[Original CarPS path migration](https://github.com/Moke110/carps/blob/ac62a177a54aa99e6ac626728af12394f3e0d940/stage2_CarPS-VAE/README.md#deploy-the-carps-stage-2-path)
