# Stage 2 shared utilities

Add `stage2/util/output/api` to the Python import path. The public
`stage2_util` package provides source-specific drug identities and representations,
indexed gene identities, seed selection, generic metrics, CSV rendering and derivation utilities.

| Module | Entry point |
| --- | --- |
| `identity` | `load_gene_table()` verifies and reads the accepted 18,164 rows in index order. `GeneIdentity` resolves global and source-specific identifiers. |
| `representations` | `RepresentationTable(source="tenv")` or `RepresentationTable(source="tahoe_bulk")` verifies the source manifest, reads the parent identity table, and selects vectors by ordered parent InChIKeys. `ordered_parents(observations)` and `verified_parent_keys(parents)` validate profile and accepted representation identities. |
| `ecfp4` | `open_ecfp4(path)` verifies the accepted binary fingerprint and parent order. `publish_ecfp4(source, candidate)` preserves frozen input provenance. |
| `seeds` | `read_seeds()` returns the accepted seed tuple. Pass a count to select a prefix. |
| `metrics` | `vector_metrics(truth, prediction)` scores paired finite vectors, including the accepted constant-vector policy. |
| `plotting`, `palette` | `read_rows(path)`, `save_figure(figure, path, **kwargs)` and the drawing functions below support rendering solely from plotted CSV values. |
| `publication` | `publication_lock(destination)` and `replace_release(...)` coordinate verified replacement of accepted releases. |
| `derivation` | `bulk_deg`, `disease_deg`, and `deg_rank` provide common DEG, rank, and worker interfaces. `r_env` locates R resources and `release` provides derivation manifest and candidate verification. |

## Plotting

The [palette figure](figure/palette.png) shows all 36 hue–depth coordinates
and the mark roles defined in [ADR 0093](../../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md).
Run `python stage2/util/figure/script/render_palette.py` from the repository
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
PYTHONPATH=stage2/util/output/api python -m unittest discover -s stage2/util/execution/tests
```
