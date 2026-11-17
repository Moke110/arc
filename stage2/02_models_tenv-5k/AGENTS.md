# Single-dose tenv 5k comparison

Read README.md, ../AGENTS.md and ADR 0074 before changing this step.
The public implementation is output/api/tenv5k; execution imports that code.
Inputs come through Bases public bounded readers; source bytes stay immutable.
The approved behavioral test seams are #226 Testing Decisions.
Build candidates under execution/checkpoints; publish only after complete acceptance.
Run one heavy task at a time under the repository resource cap and /usr/bin/time -v.
Preserve ignored inputs and outputs before any worktree cleanup.
