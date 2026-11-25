# Stage 2 perturbation prediction

## Selected gene axis

**Denoise-DEG 5k (`denoise-5k`)** is the Stage 2 prediction axis.
The [Step 04 matched-control comparison](04_tsg_axis/README.md#gene-axis-selection)
supports this choice based on held-out perturbation Delta prediction.
Across six axes and three cold-drug folds, Denoise-DEG 5k had the highest
mean Δ PCC (0.4066) and the lowest mean Δ RMSE (0.2596) and Δ MMD (0.1793).

Use the [ordered 5,000-gene list](02_gene_axis/02.2_denoise-deg/output/release/gene_axis_5000.csv)
and its [Step 03.7 tenv training release](03_data_train/03.7_tenv_denoise-5k/README.md).

## Local history

[ADR 0087](../docs/adr/0087-stage2-active-steps-and-history.md) defines the active
Steps 01–04 and shared utilities. The materialized checkout keeps superseded
raw and training releases, prior comparisons, retired Steps 05–09 and rollback
materials in `stage2/history`. The local `archive-manifest.json` records original
paths, archived paths and replay revisions. Normal readers, builds and
comparisons use the active releases. In an archive-isolated checkout with the
active inputs materialized, run `python stage2/04_tsg_axis/execution/verify_active_inputs.py`
to verify the Step 01–03 readers, projections, and fold population.

The [cutover smoke receipt](04_tsg_axis/log/323-archive-isolated-smoke.json)
records the archive-isolated acceptance.

### Deploy the ADR 0087 cutover

From the clean materialized `main` checkout, use the reviewed `dev` helper
from its worktree. Preparation checks the target Git ignore rules and detects
occupied destinations before the merge. Keep external inputs and other
Stage 2 jobs idle while installing and verifying.

```bash
STATE="$PWD/tmp/stage2-history-cutover"
python "$DEV_WORKTREE/script/stage2_history.py" prepare --root "$PWD" --state "$STATE" --target dev &&
git merge --no-overwrite-ignore --no-ff --no-edit dev &&
python script/stage2_history.py install --root "$PWD" --state "$STATE" &&
python script/stage2_history.py verify --root "$PWD" --state "$STATE"
```
