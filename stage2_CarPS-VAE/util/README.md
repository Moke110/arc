# Stage 2 shared utilities

Add `stage2_CarPS-VAE/util/output/api` to the Python import path. The public
`stage2_util` package provides source-specific drug identities and representations,
indexed gene identities, seed selection, generic metrics, CSV rendering and derivation utilities.

| Module | Entry point |
| --- | --- |
| `identity` | `load_gene_table()` verifies and reads the accepted 18,164 rows in index order. `GeneIdentity` resolves global and source-specific identifiers. |
| `representations` | `RepresentationTable(source="tenv")` or `RepresentationTable(source="tahoe_bulk")` verifies the source manifest, reads the parent identity table, and selects vectors by ordered parent InChIKeys. `ordered_parents(observations)` and `verified_parent_keys(parents)` validate profile and accepted representation identities. |
| `ecfp4` | `open_ecfp4(path)` verifies the accepted binary fingerprint and parent order. `publish_ecfp4(source, candidate)` preserves frozen input provenance. |
| `seeds` | `read_seeds()` returns the accepted seed tuple. Pass a count to select a prefix. |
| `metrics` | `vector_metrics(truth, prediction)` scores paired finite vectors, including the accepted constant-vector policy. `signed_jaccard_scores(truth, prediction, gene_indices)` returns JAC-top200 and wJAC-2000 for one paired Delta vector. |
| `comparison` | `write_paired_metrics(data, reference)` and `render_comparison(data, output, arms=..., reference=..., title=...)` publish the three canonical comparison figures. |
| `comparison_publication` | `render_receipt(stage, step, output, renderer, dependencies=...)` replays the figures twice and records the plotted CSV, PNG and renderer hashes. |
| `plotting`, `palette` | `read_rows(path)`, `save_figure(figure, path, **kwargs)` and the drawing functions below support rendering solely from plotted CSV values. |
| `publication` | `publication_lock(destination)` and `replace_release(...)` coordinate verified replacement of accepted releases. |
| `derivation` | `bulk_deg`, `disease_deg`, and `deg_rank` provide common DEG, rank, and worker interfaces. `r_env` locates R resources and `release` provides derivation manifest and candidate verification. |

## Comparison figures

Use the four shared signed Jaccard cutoffs from [ADR 0102](../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md)
at the held-out condition boundary. Each study keeps its established fold
and seed aggregation. Provide the comparison renderer with these CSVs in
its `data` directory.

| CSV | Columns |
| --- | --- |
| `learning-curve.csv` | `arm,fold,seed,epoch,train,validation` |
| `metrics.csv` | `arm,fold,seed,metric,value` |
| `metrics-paired.csv` | `arm,fold,seed,metric,reference,arm_value,reference_value,value` |

Use `comparison.METRIC_NAMES` for metric IDs. Keep undefined metric values
blank. `write_paired_metrics` derives the paired CSV from matching reference
rows. Give the renderer an ordered mapping from arm IDs to `(label, hue, depth)`
and the reference arm ID. It checks fit coverage, epoch order and paired
values before producing the three PNGs.

## Plotting

The [palette figure](figure/palette.png) shows all 36 hue–depth coordinates
and the mark roles defined in [ADR 0093](../../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md).
Run `python stage2_CarPS-VAE/util/figure/script/render_palette.py` from the repository
root to reproduce it. Use `--output PATH` to render to a different file.

- `palette.color_pair(hue, depth=2)` returns the `main` and `secondary`
  RGBA colors. Hue is in degrees and depth runs from 0 (darkest) to 3
  (lightest). The palette figure labels those depths 1 to 4.
- `plotting.bar(axis, x, height, hue=..., depth=...)` accepts one palette
  coordinate for the series or per-bar hue and depth arrays.
  `grouped_bar(axis, x, groups, width=.8)` centers each group at its x
  value and divides the total width among its series. Each series provides
  `(heights, hue, depth, label)`.
- `line`, `scatter`, `errorbar`, `range_band`, `bubble`, and `histogram`
  take one coordinate per series. Pass measurement errors as `yerr` or
  `xerr` to `errorbar` and the existing bin edges as `bins` to `histogram`.
- `plotting.replicate_dots(axis, x, values, hue=..., depth=..., spread=...)`
  draws one dot per defined replicate at the group center. Its symmetric
  spread separates dots without using fold or seed as a position variable.
- `heatmap(axis, values, hue=..., vmin=..., vmax=..., colorbar_label=...)`
  uses the four ordered depths of one hue and adds a numeric color scale.

- `plotting.replicate_bar(axis, x, values, hue=..., depth=...)` draws the
  mean, one dot per defined replicate and sample SD for at least two
  defined measurements. Use `replicate_errorbar` for the same dots and
  error bars on an explicitly requested non-bar chart.
- `plotting.learning_curve(axis, runs)` draws one training and one
  validation line per run. Pass its observed `(epochs, train, validation)`
  values and call once per model subplot.
- `plotting.focus_y(axis, groups, reference=...)` sets a padded y range
  covering finite observations and sample SD. Pass one iterable per group
  and provide a reference value when the panel needs it in view. Use the
  same collection of groups for panels that share a scale.

Use one color axis per grouping variable. Map the first variable to hue
and the second to depth when present. Replicates inherit the color of
their group. See the [Stage 2 contract](../AGENTS.md#step-contract) for
publication layout and chart defaults. Study renderers own the layout
and labels. Drawing functions return Matplotlib artists.

The default readers open the accepted gene and feature products under util.
For a separate materialized checkout, pass the directory explicitly to `load_gene_table` or `RepresentationTable`. The
feature reader checks the shipped source-specific manifest pin unless an
explicit expected manifest digest is provided. The manifest binds the frozen
input generation, extraction provenance, parent structures and vector files.
Studies retain fold coverage verification, aggregation and report tables.

From the repository root:

```sh
PYTHONPATH=stage2_CarPS-VAE/util/output/api python -m unittest discover -s stage2_CarPS-VAE/util/execution/tests
```


## TranSiGen training and evaluation

Steps 07–10 use the same runtime. `util/training/alpha_model.py` implements the fixed
CarPS-VAE alpha. `util/training/transigen.py` owns shared epoch selection,
minimum-validation checkpointing and exact replay. Compare the study-owned rows,
dose transforms, evaluation seeds and metric definitions when using these APIs.
`util/training/registered.py` handles locked phase receipts for Steps 08–10.
`util/eval/axis_delta.py` adds the weighted signed Jaccard to the shared scoring
kernels for the Step 09 and 10 full-axis comparisons. Step 09 retains its
tracked historical common-axis evidence and the accepted full-axis report.
Step 10 scores the held-out drug-and-dose intersection with its own fold calibration.


`training/transigen.py` owns registration, fit, checkpoint validation, prediction
and independent replay for Steps 04–06. Each study supplies its roster, model
builder and objective. New registrations explicitly identify architecture,
objective and weights, chemical representation, fitted dose transform and
conditioning, ordered axis digest, split, fold and training seed. The native
study records the author MSE and KL coefficients. The numbered run commands
retain their arguments and continue to require separate processes for each phase.

`eval/fixed_prediction.py` checks ordered held-out conditions, paired predictions
and their signed Delta arrays before computing condition and fold metrics.
Step 04 adds its disease-direction score and each study owns its aggregation,
training calibration, reused coordinates and publication.

The [verification record](training/verification.md) lists formal products, code
revisions and parity checks. Historical registrations retain the original
recorded code hashes. Replay them in the pinned Git checkout with the verified
input locations; new registrations bind the shared runtime.

Run one test group at a time from the repository root, for example

```bash
python -m unittest discover -s stage2_CarPS-VAE/util/training/tests/06.4_round4
python -m unittest discover -s stage2_CarPS-VAE/util/eval/tests/04_tsg_axis
```
