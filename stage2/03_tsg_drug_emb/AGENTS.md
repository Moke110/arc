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

## Fit and replay

Run `execution/smoke.py` with the accepted tenv input and representation
freeze, the pinned native TranSiGen source, a fresh destination, and
`--device cuda`. Use the transfer step's `execution/capped_phase.sh` to
collect resource evidence. The [smoke receipt](log/270-bounded-smoke.json)
records seven two-epoch origins. Preserve the separate
[resource receipt](log/270-bounded-smoke-resources.json) and raw phase evidence.

Run `execution/formal.py run` with the same input, representation and source
paths. Pass the smoke acceptance JSON and a fresh destination. The driver
runs capped phases singly, records a validation selection for each complete
fold before outer scoring, and recomputes Mean_cell on the same partitions.
Reissue the same command to resume. Use `execution/formal.py verify` with
those arguments to independently verify the completed matrix and evidence.

## Report preparation

After `execution/formal.py verify` passes, run `execution/publish_report.py`
with the formal matrix, frozen tenv inputs and representation candidate as
`--matrix`, `--inputs` and `--representations`. Give it a fresh `--destination`
and wrap the command with the same capped launcher. The report builder checks
the formal receipts and held-out identities, writes the plotted CSVs and
verifies byte-exact rendering in two fresh processes. Review the three actual
PNGs and run the offline path, tracking and integrity gates before accepting
the report.

Inspect the actual candidate PNGs at readable scale before publication. Run
`execution/report_publication.py verify-candidate` and then `publish` in the
report ticket worktree with `--candidate`, `--inputs`, `--representations` and
`--matrix`. Pass the observed figure review to `publish --visual-review`.
Publication seals the report inventory, verifies the installed files through
fresh CSV-only rendering and preserves existing records in `log/`. Run the
repository offline health gates and `verify-published` with the final payload
locations before deployment.

Install the validated feature freeze under `output/representations/`. Keep the
full registered formal matrix, checkpoints and replays together under
`execution/checkpoints/formal/` so the report receipt can resolve its original
acceptance digest. These are checkout-local payloads. The report's README,
plotted CSVs, three PNGs, acceptance receipt and run inventory are tracked.
