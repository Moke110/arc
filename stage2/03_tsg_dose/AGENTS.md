# TranSiGen dose-input optimization — round1

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md) and
[ADR0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md) before changes.

**Active execution.** Use `execution/axis_local/formal.py` and
`output/api/tsg_round1` for the 84-fit axis-local comparison. Keep Tahoe and
tenv independent. Lock cold-drug selections on validation results before
outer testing. Cold-drug-dose compares fixed methods. Run formal fits singly
under `execution/capped.sh` with wall time and 10-second RAM/GPU/WSL monitoring.

**Report publication.** Follow the ADR0078 release procedure in README.md.
Reproduce all three active PNGs from the CSVs in `figure/data` using
`figure/script/render_axis_local.py`. Complete exact replay, scientific and
visual acceptance, and the offline health checks before installation.

**Historical reconstruction.** Follow README.md and ADR0075/0076 for the
126-fit and 147-fit records. Use `execution/historical/report_126/` for the
sealed 126-fit CSV report and `execution/historical/render_v1.py` for the
12-PNG generation. Preserve the original input, model, fit and receipt bytes.
Use the approved bounded Bases readers. Accepted Stage 2 02 and Stage 3
keep independent scientific and weight contracts.
