# TranSiGen drug representation comparison

Read [the report](README.md), [Stage 2 routing](../AGENTS.md),
[ADR 0079](../../docs/adr/0079-stage2-drug-representation-step.md),
and the [execution environment](../../docs/agents/execution-environment.md).
Reuse the accepted tenv cold-drug folds and dose conditioning from the
[dose comparison](../04_tsg_dose/AGENTS.md).

## Representation input freeze

Run `execution/freeze_representations.py --help` to see the required input
paths. Use the accepted Round 1 tenv freeze at
`../04_tsg_dose/execution/checkpoints/round1-axis-local-inputs/tenv`,
the dose `output/api`, the vendored KPGT source and its frozen base checkpoint,
and the official UniMol v1 all-hydrogen checkpoint and molecular dictionary.
Run the extraction under the dose step's `execution/capped.sh` with
`/usr/bin/time -v`. Use `e3fp==1.2.7` and `unimol_tools==0.1.5` in the
primary Python environment and pass the dedicated KPGT environment's Python.
For example, from the repository root, with those interpreters and UniMol
weights configured:

```sh
stage2/04_tsg_dose/execution/capped.sh /usr/bin/time -v \
  "$ARC_PYTHON" stage2/03_tsg_drug_emb/execution/freeze_representations.py \
  --inputs stage2/04_tsg_dose/execution/checkpoints/round1-axis-local-inputs/tenv \
  --dose-api stage2/04_tsg_dose/output/api \
  --kpgt-python "$KPGT_PYTHON" \
  --kpgt-source drug_emb_model/KPGT/source \
  --kpgt-weights drug_emb_model/KPGT/source/models/pretrained/base/base.pth \
  --unimol-weights "$UNIMOL_WEIGHTS/mol_pre_all_h_220816.pt" \
  --unimol-dict "$UNIMOL_WEIGHTS/mol.dict.txt" \
  --output "$REPRESENTATION_CANDIDATE"
```

The candidate contains one structure ledger, seven native-scale float32
matrices, the per-parent E3FP conformer counts, and a digest manifest. The
extractor writes `coverage_gaps.json` instead of a manifest if it cannot
cover the complete frozen roster. Verify through
`tsg_drug_emb.representations.RepresentationInputs` before passing the
candidate to training.
