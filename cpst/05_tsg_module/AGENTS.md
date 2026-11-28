# Step 05 TranSiGen module workflow

Read the [study report](README.md), [accepted design](../../docs/adr/0110-cpst-transigen-module-rerun.md)
and [CPST instructions](../AGENTS.md). Work from the repository root in
`arc`. Follow the [execution limits](../../docs/agents/execution-environment.md)
and sample WSL, Windows and GPU every ten seconds.

Use `execution/matrix_measured.py` inside a 12 GiB RAM, 8 GiB swap systemd
scope, with CPU thread counts capped at two. It registers all 125 coordinates
before fitting and resumes from phase receipts. Pass `--register-only`
to stop after verifying the complete registration ledger. Set `--training`,
`--representations` and `--source` to the materialized Denoise-DEG 5k,
KPGT and native TranSiGen sources. Place the matrix, separate audit logs,
temporary directory and score report under a fresh directory in `tmp/`.

Use the matrix with `--smoke` for the bounded gate. It uses the complete
membership of the first fold, seed 17, and two epochs for each module.
Run `execution/audit_measured.py objectives --smoke` against the matrix,
and `execution/verify_checkpoint.py` in a capped CUDA scope for an independent
selected-checkpoint loss and validation-prediction replay.

Run `execution/evaluate_measured.py score` and `verify` in separate processes
with the matrix, training release, HVG-DEG 5k bandwidth release, report output
and scratch paths. Pass `--smoke` only for the bounded matrix. These commands
derive held-out scores from sealed predictions and matched control and verify
the selected checkpoint receipts. Run `execution/audit_measured_scores.py`
on the report for an independent audit of the coordinate-paired tables.

Build figures with `execution/build_comparison.py --report <report> --output <directory>`
using an unused output directory, with `--smoke` for the bounded report. The
receipt records two independent CSV-only renders. Use `figure/script/render.py`
for CSV-only replay with `--data` and `--output`.

## Accept and replay the installed publication

Score the complete formal matrix first. Then run
[`publish_module.py`](execution/publish_module.py) with `build`,
`--publication` set to an unused candidate directory in `tmp/`, and the
sealed matrix, matrix audit, and score report paths. Supply
`--coordinator-samples` with the ten-second resource log and
`--code-revision` with the revision used to run the matrix. The command
checks the preregistration, repeats test scoring in a fresh process, audits
the objective traces and coordinate-matched tables, and renders all three
figures twice from CSV. Its formal receipt seals the accepted outputs and
resource evidence.

Install the accepted matrix, matrix audit and score report as
`output/formal/`, `output/formal-audit/` and `output/report/` respectively.
The score report retains the condition-level results and fits retain the
checkpoint, test predictions and separate-process replay evidence. Publish
the projected tables, figures and receipts at the step's tracked paths.

In the materialized checkout, rerun verification with the originally pinned
inputs and environment. Set `TMPDIR` and `MPLCONFIGDIR` to directories under
checkout-root `tmp/`, and run under the execution guide's memory limits.

```bash
step="$PWD/cpst/05_tsg_module"
work=$(mktemp -d "$PWD/tmp/step05-installed.XXXXXX")
mkdir "$work/matplotlib"
export TMPDIR="$work" MPLCONFIGDIR="$work/matplotlib"
python "$step/execution/publish_module.py" verify \
  --matrix "$step/output/formal" \
  --matrix-audit "$step/output/formal-audit" \
  --report "$step/output/report" --scratch "$work"
```
