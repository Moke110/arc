# TranSiGen drug representation comparison

Read [the report](README.md), [Stage 2 routing](../AGENTS.md),
[ADR 0079](../../docs/adr/0079-stage2-drug-representation-step.md),
and the [execution environment](../../docs/agents/execution-environment.md).
Reuse the accepted source-local cold-drug folds and dose conditioning from the
[dose comparison](../04_tsg_dose/AGENTS.md).

## Raw-first comparison

Use `execution/raw_first.py` for new comparison preparation, fixed-prediction
evaluation and CSV-only plotting. The input is the source directory emitted by
the [dose step's raw-first preparation](../04_tsg_dose/README.md#raw-first-inputs-evaluation-and-figures).
Follow the [active commands](README.md#raw-first-inputs-evaluation-and-figures).
The instructions below reproduce the sealed historical comparisons.

## Representation input freeze

Run `execution/freeze_representations.py --help` to see the required input
paths. Use the accepted Round 1 source directory under
`../04_tsg_dose/execution/checkpoints/round1-axis-local-inputs/`.
Use `tenv` for the 459-parent roster or `tahoe_bulk` for the 352-parent
roster. Use the dose `output/api`, the vendored KPGT source and its frozen base checkpoint,
and the official UniMol v1 all-hydrogen checkpoint and molecular dictionary.
Run the extraction under the dose step's `execution/capped.sh` with
`/usr/bin/time -v`. Use `e3fp==1.2.7` and `unimol_tools==0.1.5` in the
primary Python environment and pass the dedicated KPGT environment's Python.
Get the two UniMol v1 assets from [DP Technology’s official upload](https://huggingface.co/dptech/Uni-Mol-Models/commit/b6904ae06b1257be231c6fd44795de8d1cde03ec).
Set `UNIMOL_WEIGHTS` to a local directory containing both files. To download
them into a new directory, run the following with the Hugging Face CLI.

```sh
hf download dptech/Uni-Mol-Models mol_pre_all_h_220816.pt mol.dict.txt \
  --revision b6904ae06b1257be231c6fd44795de8d1cde03ec \
  --local-dir "$UNIMOL_WEIGHTS"
```

From the repository root, verify both assets against the tracked freeze
receipt before running the extractor. The same check applies to existing
files supplied by `unimol_tools`.

```sh
(
  set -euo pipefail
  receipt="$PWD/stage2/03_tsg_drug_emb/log/representation-freeze.json"
  cd -- "$UNIMOL_WEIGHTS"
  jq -er '.provenance.UniMol |
    (.checkpoint, .dictionary) | "\(.sha256)  \(.file)"' "$receipt" | sha256sum --check --strict
)
```

Continue after both SHA-256 checks pass. From the repository root, with the
two interpreters configured, set `INPUT_SOURCE` to `tenv` or `tahoe_bulk`:

```sh
stage2/04_tsg_dose/execution/capped.sh /usr/bin/time -v \
  "$ARC_PYTHON" stage2/03_tsg_drug_emb/execution/freeze_representations.py \
  --inputs "stage2/04_tsg_dose/execution/checkpoints/round1-axis-local-inputs/$INPUT_SOURCE" \
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

Run `execution/smoke.py` with the accepted source-local input and representation
freeze, the pinned native TranSiGen source, a fresh destination, and
`--device cuda`. Use the transfer step's `execution/capped_phase.sh` to
collect resource evidence. The [smoke receipt](log/270-bounded-smoke.json)
records seven two-epoch origins. Preserve the separate
[resource receipt](log/270-bounded-smoke-resources.json) and raw phase evidence.

Run `execution/formal.py run` with the same input, representation and source
paths. Pass the smoke acceptance JSON and a fresh destination. The driver
runs capped phases singly, records a validation selection for each complete
fold before outer scoring, and recomputes Mean_cell on the same partitions.
Reissue the same command to resume after a phase with complete telemetry
or after an attempt that left no phase directory. Use
`execution/formal.py verify` with those arguments to independently verify the
completed matrix and evidence.

### Recover an interrupted phase

A phase can write its directory before the capped launcher finishes writing
resource evidence. In that case the driver rejects reuse until the incomplete
phase is quarantined. Keep the formal plan and all accepted evidence intact.

1. Stop the driver before moving files. Confirm that its matching
   `arc-stage2-transfer-<phase-attempt>.scope` is inactive in
   `systemctl --user list-units --type=scope` and that no process can still
   write the phase directory or telemetry.
2. Use the phase identifier in the failure to locate the exact directory
   and corresponding attempt directories under `resources/phases/`. The
   run directories are `runs/fold-<fold>/<arm>/`, with `train` mapped to
   `fit`, `replay-validation` to `validation-replay`, `evaluate` to
   `test` and `replay-test` to `test-replay`. The Mean_cell phase directory
   is `mean-cell/fold-<fold>`.
3. Quarantine only if the destination has no `acceptance.json`, the phase
   has no successful attempt ending in `status=0` and no later phase has
   a receipt. The driver processes folds in order. Within a fold it runs
   each arm’s fit and validation replay, locks the fold selection, runs
   each arm’s test and test replay, then recomputes Mean_cell. A fold
   selection already binds the fits and validation replays. If any
   evidence binds the interrupted phase, preserve the entire destination
   and start a new one.
4. Hold an exclusive lock on the destination’s `.orchestration.lock`
   while moving just the interrupted phase directory and its failed or
   incomplete telemetry attempts to a fresh quarantine directory beside
   the destination. Keep their original relative paths and release the
   lock before retrying. For example, use `flock -n "$MATRIX/.orchestration.lock" bash`
   to open a locked shell for the checked moves, then exit that shell.
5. Reissue the original `formal.py run` command with the same destination
   and inputs. A new capped attempt recreates that phase. Run `formal.py
   verify` after the matrix completes.

## Report preparation

After `execution/formal.py verify` passes, run `execution/publish_report.py`
with the formal matrix, inputs from the same source and representation candidate as
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

For Tahoe, use a separate candidate during coverage and fit. Stop before any
fit unless all seven matrices pass on the complete 352-parent roster. Install
the accepted local freeze at `output/representations/tahoe/` and the verified
fit payload at `execution/checkpoints/formal-tahoe/`. The Tahoe report places
its receipt and results in `log/tahoe/`, its plotted tables in
`figure/data/tahoe/` and its three prefixed PNGs in `figure/`. Render using
`python figure/script/render.py --source tahoe_bulk --prefixed --output <new-directory>`
from this step.

The original tenv runtime and report renderer are pinned by the historical
receipts. Check out commit `107bda7` in a separate worktree and pass the
existing absolute input, representation and run paths to reproduce that
byte version independently.

Publishing the paired report renames the three original `figure/<name>.png`
files to `figure/tenv_<name>.png` and moves the sealed `README.md` to
`log/tenv/README.md` byte for byte. The tenv receipt remains at
`log/acceptance.json`; `verify-published` resolves those historical names via
the migration mapping and checks the same digests and CSV-only replay.
