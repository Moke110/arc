# TranSiGen Tahoe-to-tenv transfer

[Specification #255](https://github.com/Moke110/arc/issues/255) evaluates
Tahoe pretraining and tenv fine-tuning across three tenv cold-drug-dose
folds. Formal acceptance and the paired comparison are pending.

The study uses the Step 03 32-dimensional Latent embedding and axis-local
`G_5000` expression. Tahoe keeps replicates of each parent, cell line and
dose condition together and holds out each eligible condition in at most
one validation fold. Each validation parent, cell line and dose has
independent training coverage. Three Tahoe fits restore the strict
minimum treated-expression validation MSE. The lowest fold MSE selects
one immutable checkpoint, with equal MSE broken by the earlier fold ID.

## Tahoe reproduction

Use the repository's `arc` environment. If the Step 03 input freeze is
unavailable, rebuild it through Step 03's documented preparation command.
Set `PAIR` to its output directory, `TAHOE` to an empty candidate path
and `SOURCE` to the pinned native TranSiGen source checkout.

```bash
STEP=stage2/04_tsg_transfer
CAP=stage2/03_tsg_dose/execution/capped.sh
"$CAP" /usr/bin/time -v python "$STEP/execution/tahoe.py" prepare \
  --inputs "$PAIR/tahoe_bulk" --candidate "$TAHOE"
for fold in 0 1 2; do
  for action in fit replay; do
    "$CAP" /usr/bin/time -v python "$STEP/execution/tahoe.py" "$action" \
      --inputs "$PAIR/tahoe_bulk" --candidate "$TAHOE" \
      --native-source "$SOURCE" --fold "$fold" --device cuda
  done
done
"$CAP" /usr/bin/time -v python "$STEP/execution/tahoe.py" select \
  --inputs "$PAIR/tahoe_bulk" --candidate "$TAHOE"
```

Use a separate candidate for a bounded smoke. Add `--smoke-rows 128
--smoke-epochs 2` to fits and `--smoke-rows 128` to replays and selection.
The split CSV accounts for training-only profiles. Fit receipts bind the
training-only standardizer, ordered partitions, checkpoint, learning
history and exact validation replay.
