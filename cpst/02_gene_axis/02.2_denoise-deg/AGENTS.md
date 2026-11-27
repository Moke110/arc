# Denoise-DEG prediction axes

Read the [prediction-axis routing](../AGENTS.md) and the [selection report](README.md).
Downstream callers use [`prediction_axis.open_axis`](../02.1_hvg-deg/output/api/prediction_axis.py)
with `denoise-1k` through `denoise-5k`. Use `axis.project(corpus, "tenv")` for
bounded native counts in the published gene order.

## Reproduction

1. Check the three inputs against the [release manifest](output/release/manifest.json)
   before rebuilding. Use the installed [Step 01 corpus](../../01_data_raw/README.md)
   manifest and its native tenv payload, the [02.1 release manifest](../02.1_hvg-deg/output/release/manifest.json)
   with its eligible disease DEG1000 lists, and the frozen [gene identities](../../util/output/gene-identity/genes.csv).
   Resolve the historical `stage2/` manifest keys under `cpst/`.
2. For a full rebuild, copy `execution/` to a candidate step under repository `tmp/`.
   Run [diagnose.py](execution/diagnose.py) and then [diagnose_scores.py](execution/diagnose_scores.py)
   with `--repo` set to the input checkout and `--step` to the candidate. Follow the
   [execution limits](../../../docs/agents/execution-environment.md). The scoring command
   uses the [plate worker](execution/plate_deg.R) with limma-voom and native library
   depths. Use `ARC_R_HOME` to select the R installation through the [shared derivation API](../../util/README.md).
3. Run [build_axes.py](execution/build_axes.py) with the same `--repo` and `--step`.
   For an evidence-only replay, start here with a copy of `figure/data/` in the
   candidate. Keep same-plate matched DMSO controls, one tenv vote per gene and
   parent–dose after the top-DEG cap, and one disease vote per complete eligible
   DEG1000 list. Apply the exclusion union before ranking. Break score ties by
   frozen gene index and write the five nested memberships in ascending index order.
4. Replay the figures with [render.py](figure/script/render.py) using `--data` for
   the plot tables and `--output` for the candidate PNGs. The renderer loads the
   [shared plotting API](../../util/README.md#plotting).

## Acceptance

Use the [release tests](execution/tests/test_disease_evidence.py) to check nesting,
gene order and coverage, with the full eligible DEG1000 lists as denominators.
Check input, evidence and output digests, then read all five axes and project
bounded tenv rows through the public API. Use the [step receipt](log/acceptance.json)
for selection and native replay and the [rendering receipt](../log/342-figure-render.json)
for current PNG, data and renderer digests. Promote the verified `output/release/`
and matching figure evidence under the [stage publication rules](../../AGENTS.md).
