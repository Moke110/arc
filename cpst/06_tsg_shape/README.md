# TranSiGen architecture studies

The numbered commands use the [shared TranSiGen training and evaluation APIs](../util/README.md#transigen-training-and-evaluation).

[Step 06.1](06.1_round1/README.md) compares ten expression architectures
on the fixed Stage 2 gene axis and training objective. The encoder-width
extension adds six fits and reuses the 24 accepted fits. The original
architecture supplies the paired reference on each frozen fold.

[Step 06.2](06.2_round2/README.md#formal-results) compares direct Δ readouts,
drug–dose conditioning, encoder sharing and latent residuals. Seven new
architectures share the three frozen folds and the three Step 06.1
original fits supply the paired references.

[Step 06.3](06.3_round3/README.md#formal-results) compares signed Δ readouts,
normalized FiLM and their combinations on the same three cold-drug folds
across four training seeds. Its 72-coordinate report compares each
combination with its components and the original within fold and seed.

[Step 06.4](06.4_round4/README.md#formal-results) compares measured
encoder targets and signed Δ decoder sharing across the same three
cold-drug folds and three training seeds.

## Rebuild the comparison figures

Each study publishes its learning curves, held-out metrics and differences
from its established reference within fold and seed. The weighted signed
Jaccard score uses the verified held-out predictions and matched controls.
The other scores come from the accepted reports. Use the current figure
receipts to replay the PNGs from their plotted CSVs.

From the materialized checkout in the `arc` environment, regenerate a
candidate from the accepted fit and report releases. Use `--step` to
select a single comparison and `--evidence-root` to read the materialized
releases from a separate checkout.

```bash
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-comparisons.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
python cpst/06_tsg_shape/execution/build_comparison.py \
  --output "$WORK/publication"
```
