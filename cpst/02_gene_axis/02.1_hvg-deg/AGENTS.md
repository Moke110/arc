# HVG-DEG axis method

Read the [experiment summary](README.md), [selection decision](../../../docs/adr/0073-full-source-tier-first-deg1000-and-four-axes.md) and [axis/view boundary](../../../docs/adr/0083-stage2-corpus-axis-view-boundaries.md) before changing selection or projection.

## Verify and publish

Run from the checkout root in the `arc` environment. The [axis manifest](output/release/manifest.json) binds the fixed 18,164-gene identity, bundled tenv HVG ranking, and tier-first disease DEG1000 evidence. Verification recomputes the nested 2k, 3k, 4k and 5k axes and checks the accepted gene order, symbols and hashes.

```bash
STEP=cpst/02_gene_axis/02.1_hvg-deg
CORPUS=cpst/01_data_raw/output
python "$STEP/execution/build_axes.py" verify
python -m unittest discover -s "$STEP/execution/tests"
python "$STEP/execution/build_axes.py" verify --corpus "$CORPUS"
```

Use the first command for axis-only verification on a fresh checkout. For release acceptance, materialize the [native corpus](../../01_data_raw/AGENTS.md) and run the last command. It checks complete common-gene coverage for `tenv`, `tahoe_bulk` and every other source declared complete. Follow the [WSL resource limits](../../../docs/agents/execution-environment.md) when opening the corpus.

To rehearse publication, give [`build_axes.py`](execution/build_axes.py) the `publish` action with `--corpus`, `--source-root` pointing to the checkout with the accepted axis evidence, and `--output` under the checkout’s `tmp/`. Keep historical provenance byte for byte. Publication owns axis and evidence tables. The [training data step](../../03_data_train/AGENTS.md) owns derived expression views.

## Projection and figure replay

Open [the axis API](output/api/prediction_axis.py) with `prediction_axis.open_axes().axis(size)` and project through `axis.project(corpus, source_name)` while the native corpus is open. Use `counts(start, stop)` or `batches(batch_size)` to keep reads bounded and preserve integer counts in prediction-axis order.

For the coverage figure, run [`render_coverage.py`](figure/script/render_coverage.py) with `--data "$STEP/figure/data"` and `--output` under the checkout’s `tmp/`. Compare the PNG hash with the [figure receipt](../log/342-figure-render.json). The renderer uses the [shared plotting API](../../util/README.md#plotting).
