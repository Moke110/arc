# CarPS-Transcriptome prediction axes

Read the [step evidence](README.md) and the [CarPS-Transcriptome release rules](../AGENTS.md). Step 02 publishes nested, ordered subsets of the fixed 18,164-gene common axis and the evidence used to select them. [Step 01](../01_data_raw/AGENTS.md) owns native counts. [Step 03](../03_data_train/AGENTS.md) owns pooling, normalization and training views. [ADR 0083](../../docs/adr/0083-stage2-corpus-axis-view-boundaries.md) · [ADR 0085](../../docs/adr/0085-stage2-raw-and-training-data.md)

## Select and verify axes

- **HVG-DEG.** Verify the four released axes from the frozen gene identity, bundled tenv HVG ranking and eligible disease DEG1000 lists. For acceptance or publication, use the materialized Step 01 corpus and check complete common-axis coverage for tenv, Tahoe bulk and every other source declared complete. The [02.1 guide](02.1_hvg-deg/AGENTS.md#verify-and-publish) owns the commands and candidate publication. [HVG-DEG manifest](02.1_hvg-deg/output/release/manifest.json)
- **Denoise-DEG.** Use the 02.1 disease evidence, frozen gene identity and pinned Step 01 corpus manifest. The evidence-only replay starts with accepted diagnostics and uses the corpus manifest. The full rebuild also reads the native tenv counts and runs plate-wise limma-voom against matched DMSO. Preserve the exclusion union, vote rules, tie breaks and ascending gene-index order. Use the [02.2 guide](02.2_denoise-deg/AGENTS.md#reproduction) for both paths and its [acceptance checks](02.2_denoise-deg/AGENTS.md#acceptance) before promotion. [Denoise-DEG manifest](02.2_denoise-deg/output/release/manifest.json)

## Read and project axes

From the checkout root, put `cpst/02_gene_axis/02.1_hvg-deg/output/api` on `PYTHONPATH`. Use `prediction_axis.open_axis(axis_id)` with `hvg-2k` through `hvg-5k` or `denoise-1k` through `denoise-5k`. It verifies the release schema, common-axis identity, axis digest and ordered membership. Reuse the canonical ID when selecting the corresponding tenv training view or running an axis comparison. [Prediction-axis API](02.1_hvg-deg/output/api/prediction_axis.py)

For projection, bind `axis.project(corpus, source_name)` while the native corpus is open. The source must cover the complete common axis. Read bounded windows with `counts(start, stop)` or `batches(batch_size)` to preserve integer counts and prediction-axis order. Apply the [WSL execution limits](../../docs/agents/execution-environment.md) before opening the corpus.

## Verify the current figures

From the checkout root, render both methods from their published plot tables using the `arc` environment and the versions in the [render receipt](log/342-figure-render.json). Use the [shared plotting API](../util/README.md#plotting). Keep size hues and disease depths consistent across both coverage figures and show the full 0–100% range.

```bash
AXES=cpst/02_gene_axis
mkdir -p tmp
RENDER_DIR=$(mktemp -d "$PWD/tmp/step02-figures.XXXXXX")
HVG=02.1_hvg-deg/figure
DENOISE=02.2_denoise-deg/figure
python "$AXES/$HVG/script/render_coverage.py" \
  --data "$AXES/$HVG/data" --output "$RENDER_DIR/$HVG"
PYTHONPATH="cpst/util/output/api${PYTHONPATH:+:$PYTHONPATH}" \
  python "$AXES/$DENOISE/script/render.py" \
  --data "$AXES/$DENOISE/data" --output "$RENDER_DIR/$DENOISE"
sha256sum "$RENDER_DIR"/*/figure/*.png
```

Match all five PNG hashes to the receipt and check the receipt's plot-table and renderer hashes against the checkout. For the stage-wide check, use the [figure verifier](../util/figure/script/verify_report_figures.py), which verifies those hashes and redraws each figure twice in scratch.
