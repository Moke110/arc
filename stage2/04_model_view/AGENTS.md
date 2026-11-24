# Single-dose tenv 5k comparison

Use the [historical evidence index](../03_data_train/03.1_hvg-deg-5k/historical/README.md) for sealed records.

Read README.md, ../AGENTS.md and ADR 0074 before changing this step.
Use `execution/raw_first.py` and `output/api/tenv5k_raw_first` for well-level study input
preparation and evaluation. Render an explicit summary CSV with `figure/raw_first/render.py`.
The approved behavioral test seams are #226 Testing Decisions.
Build candidates under execution/checkpoints; publish only after complete acceptance.
Run one heavy task at a time under the repository resource cap and /usr/bin/time -v.
Preserve ignored inputs and outputs before any worktree cleanup.
