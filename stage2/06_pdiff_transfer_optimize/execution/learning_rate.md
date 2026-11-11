# Base-learning-rate intervention

The seed 17 cold-drug fold 0 validation intervention changes target AdamW base
learning rate from 0.0002 to 0.0001. Each pair uses the original Tahoe EMA,
bridge, recipient and native recipe. Checkpoints remain under `output/checkpoints/`.

The plan freezes common optimizer steps separately from native loss selection.
`native-selection.json` records the selected EMA; `diagnostic.json` identifies
the fixed-step EMA. A native stop before the common step continues the same
optimizer, scheduler, random stream and EMA under the original wall-time limit,
keeping diagnostic observations out of native selection history.

AP@100 score commands, score seals and promotion thresholds are retired.
Historical retrieval results remain only in Git; the new metric does not
retroactively authorize new screening seeds. Frozen predictions now use:

```bash
python stage2/evaluation/publish.py --step 6
python stage2/evaluation/evaluate.py replay --step 6
```

The preserved CPU fixture checks exact native EMA and reload parity. Its
recorded run used two threads, 10 GiB memory and 8 GiB swap, took 28.39 seconds
and peaked at 4,182,272 KiB RSS.
