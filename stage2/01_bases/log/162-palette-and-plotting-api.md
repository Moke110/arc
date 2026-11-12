# #162 — stage palette and plotting API

Ticket #162 (parent #158), branch `codex/tenv-162` off `dev` (`8e5fbcf`, which
contains the #160 evaluation core). Authority: ADR 0063 for the palette, the two
figure layouts and the undefined rule; ADR 0062 and `stage2/AGENTS.md` for the
product/figure contract. No scientific metric is recomputed here: both drawing
APIs present values that the step computed elsewhere, and they return exactly
the values they drew.

## Palette product — `output/palette/`

`palette.csv` is the two-dimensional table: one `shade` column followed by the
six hue families (`gray,blue,green,orange,purple,magenta`), nine rows from the
darkest/most saturated shade 0 to the lightest shade 8. `manifest.json` records
the generation recipe, the colour roles, the table digest and the digest of the
public module that generates it. The table is a pure function of the recipe in
`output/api/bases/palette.py` (`RECIPE`, `color`, `table`, `publish`); no file is
read to compute a colour, so a renderer needs only recorded plot values plus the
public implementation.

Recipe: per shade `s`, `L = 0.22 + 0.07*s` and `S = 0.90 - 0.06*s` (S = 0 for
gray), converted from HLS with `colorsys.hls_to_rgb` and rounded to 8-bit.
Families carry the base hues gray achromatic, blue 214°, green 137°, orange 28°,
purple 276°, magenta 318°. Golden cells: blue/2 `#1452a3`, green/2 `#14a33d`,
purple/2 `#6a14a3`, magenta/8 `#deafd0`.

Roles the two figures need (`palette.ROLES`, recorded in the manifest): `text`,
`grid`, `undefined_marker`, `undefined_label` (fixed gray cells); `line`,
`bar_edge`, `dot` (entry hue, shade 2); `band`, `bar_fill` (entry hue, shade 6);
`train` (blue family) and `validation` (green family) at the run shade. Entry
order without an explicit hue is `blue, orange, green, purple, magenta, gray`, so
a deterministic reference takes gray last; four runs use shades `2, 4, 6, 8`.

The product files are force-added: `.gitignore` excludes
`stage2/01_bases/output/*` and re-includes only `output/api/**`; the established
#156 practice is to force-add compact evidence. No `.gitignore` rule was changed.

## Plotting API — `output/api/bases/plots.py`

Added alongside the untouched retired three-AUPRC `render_retrieval`/
`replay_retrieval` panel, which still reproduces its committed PNG byte for byte
(sha256 `35d14266…`, re-checked here):

- `metrics_values(entries)` → the complete eight-panel plotted values;
  `draw_metrics(values, axes)` draws them and returns them; `render_metrics(...)`
  saves the PNG; `replay_metrics(csv, png)` redraws from the CSV alone.
- `learning_curve_values(models)` → the complete learning-curve values;
  `draw_learning_curve`, `render_learning_curve`, `replay_learning_curve`.
- `plot_rows(values)` / `write_plot_values(values, path)` /
  `read_plot_values(path)` are the single long-form recorded-value schema shared
  by both figures (`record` = figure/panel/series/point; `role` names what was
  drawn: `bar_mean`, `bar_run`, `curve_mean`, `curve_min`, `curve_max`,
  `curve_run`, `line`). Colours are not stored as hex: the recorded hue, shade
  and role resolve through `bases.palette`, so a colour cannot drift from the
  stage palette.

`metrics.png` values: two rows × four columns in the order `PCC, Spearman, RMSE,
R2, Precision@k, Spearman@k, max_k@P80S40, Recall@max_k`. `Precision@k` and
`Spearman@k` draw `k` with the cross-run mean line and a min–max band (counts of
defined runs recorded per `k`); the six direct/threshold panels draw a bar of the
cross-run mean with one dot per run. An entry with no defined run for a metric is
NaN with a reason and is drawn as an `×` plus an `undefined` label — never as a
zero-height bar or a zero line.

`learning_curve.png` values: one subplot per learnable model, epoch on the
horizontal axis and loss on the vertical. Training is the blue family,
validation the green family, and one shade identifies one `(fold, seed)` run in
both. A run with no recorded trajectory is recorded as undefined and annotated in
the panel; it is never drawn as an empty-at-zero curve.

## Fixtures (synthetic values)

Committed under `execution/tests/fixtures/` because `figure/` belongs to the
step's own tickets; these are test evidence, not report figures.

| file | sha256 |
|---|---|
| `metrics_fixture.png` | `aa00c92c92d9ff3e3b4fb98d727294097ac78a7fd8efdde9352a06719130595e` |
| `metrics_fixture.csv` | `071c34b4b6a93bebebe8fe6a1a18da281c2ac792422d68642a8479b15760d4ee` |
| `learning_curve_fixture.png` | `62127d02e1f54ce29a7fe755a235067c7d0cf2098278d37fca33a1468caf63c0` |
| `learning_curve_fixture.csv` | `58bd73c7540ed0060a7f0dc36f823a9d85da0ee4304f5277f1b8435bd3889110` |

Regenerate with (from `stage2/01_bases/execution`):

```bash
MPLCONFIGDIR=/tmp/mplcache /home/chang/miniconda3/bin/python -c "
import sys; sys.path[:0]=['tests','.']
import public_api
from pathlib import Path
from test_metrics_panels import synthetic_entries
from test_learning_curve_panels import synthetic_models
from bases import plots
fx = Path('tests/fixtures')
plots.write_plot_values(plots.metrics_values(synthetic_entries()), fx/'metrics_fixture.csv')
plots.render_metrics(synthetic_entries(), fx/'metrics_fixture.png')
plots.write_plot_values(plots.learning_curve_values(synthetic_models()), fx/'learning_curve_fixture.csv')
plots.render_learning_curve(synthetic_models(), fx/'learning_curve_fixture.png')"
```

## Acceptance evidence (#162 checkboxes)

| checkbox | evidence |
|---|---|
| Palette is a 6×9 hex table plus recipe covering every role | `output/palette/palette.csv`, `manifest.json`; `test_palette.py` (table shape, monotone lightness, hue bounds, golden cells, recipe, all 11 roles resolve) and the publish-reproduces-committed-bytes test |
| Indirect panels mean line + min–max band, other panels mean bars with one dot per run | `test_metrics_panels.py::test_indirect_panels_are_cross_run_mean_lines_with_min_max_bands` and `::test_configured_and_threshold_panels_are_bars_of_the_cross_run_mean`; artist-level checks in `MetricsDrawingTests` |
| Entry hue colours line/edge/dots, lighter shade fills band/bar | `test_metrics_panels.py::test_bar_panels_use_the_entry_hue_edge_and_lighter_fill_with_one_dot_per_run` and `::test_indirect_panels_draw_a_band_a_mean_line_and_explicit_undefined_marks` |
| Learning curve blue/green with one shade per `(fold, seed)` | `test_learning_curve_panels.py::test_train_blue_and_validation_green_share_one_shade_per_run` and `::test_each_run_draws_blue_training_and_green_validation_at_its_shade` |
| Undefined drawn explicitly, never zero | `test_metrics_panels.py::test_undefined_bar_is_nan_with_a_reason_and_never_zero` and `::test_undefined_entry_is_annotated_and_drawn_without_a_zero_bar`; `test_learning_curve_panels.py::test_undefined_run_is_recorded_as_undefined_with_a_reason` |
| Fixture figure from synthetic values whose returned values equal what was drawn | fixture CSV/PNG above; `test_metrics_panels.py::test_draw_returns_exactly_the_values_it_drew`, `::test_recorded_values_round_trip_and_replay_byte_identical`, `::test_committed_fixture_is_the_synthetic_values_and_replays_byte_identical` |
| CSV-only reproduction needs no step `output/` | `replay_metrics`/`replay_learning_curve` read only the recorded CSV and the public palette/drawing code; the byte-identity tests prove presentation is a function of recorded values |

## Deviations and follow-ups

- Issue #162's summary sentence groups `Precision@k`, `Spearman@k` and the
  `max_k`/`Recall` panels under "plot `k`". Implemented per ADR 0063 instead:
  only `Precision@k` and `Spearman@k` are curves over `k`; `max_k@P80S40` and
  `Recall@max_k` are per-run scalars and are bars (ADR: "the six direct and
  threshold panels are bar charts"). `precision`/`spearman` are the two arrays in
  `aggregate_indirect`; `max_k`/`recall` are scalars there by construction.
- `output/palette/` is force-added (see above); no ignore rule changed.
- The step `README.md`, `AGENTS.md` and `docs/data-inventory.md` still do not
  list `output/palette/`; those files are outside this ticket's write scope and
  should be updated when the step's own figure ticket lands.
- The fixture PNGs live under `execution/tests/fixtures/`, not `figure/`,
  because `#162` owns no report figure; #163's renderer will consume this API for
  `figure/metrics.png` and `figure/learning_curve.png`.
- PNG byte identity is asserted only within one environment (matplotlib 3.10.9,
  Agg, dpi 160), the same standard the committed diagnostic figure uses.
