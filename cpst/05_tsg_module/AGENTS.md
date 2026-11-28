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
