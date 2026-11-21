# TranSiGen dose-input optimization — round1

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md) and
[ADR0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md) before changes.

**tenv KPGT successor.** Follow [ADR0080](../../docs/adr/0080-tenv-kpgt-dose-round1.md).
Prepare with `execution/kpgt/prepare.py` from the sealed Step 04 tenv input
freeze and Step 03 representations. Verify with `execution/kpgt/smoke.py`
under the repository caps. Run `execution/kpgt/formal.py` and `baseline.py`
singly with 10-second RAM/GPU/WSL monitoring, then accept via `publish.py`.
Publish its figures in `figure/` and verify them against the KPGT acceptance receipt.

**ADR0078 execution.** Use `execution/axis_local/formal.py` and
`output/api/tsg_round1` for the 84-fit axis-local comparison. Keep Tahoe and
tenv independent. Lock cold-drug selections on validation results before
outer testing. Cold-drug-dose compares fixed methods. Run formal fits singly
under `execution/capped.sh` with wall time and 10-second RAM/GPU/WSL monitoring.

**ADR0078 baselines and report.** Use `execution/axis_local/baselines.py` for
the 24 baseline slots and `publish_baselines.py` for successor acceptance.
Preserve the accepted gray presentation in
`execution/historical/report_baselines_gray/`. Complete fresh-process replay, scientific and visual
acceptance and offline checks before deploying. For sealed 84-fit and baseline
report verification, use the archived inventories described in [README.md](README.md).

**Historical reconstruction.** Follow README.md and ADR0075/0076 for the
126-fit and 147-fit records. Use `execution/historical/report_126/` for the
sealed 126-fit CSV report and `execution/historical/render_v1.py` for the
12-PNG generation. Preserve the original input, model, fit and receipt bytes.
Use bounded Bases readers. Step 02 and Stage 3 retain independent contracts.
