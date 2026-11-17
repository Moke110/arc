# All-dose TranSiGen dose-input comparison

Read README.md, ../AGENTS.md and ADR0075 before changes. Public implementation
is output/api/tsg_dose; execution imports it. Read inputs through bounded Bases
readers. Source bytes and accepted02 products stay unchanged.

Tests use #232's approved input-freeze, learner/predict/replay and score/accept
boundaries. Formal fits run singly under capped.sh and /usr/bin/time -v with
available-memory preflight and10-second RAM/GPU/WSL monitoring. Keep validation
and test phases separate: lock all-candidate validation choice before outer test.
Candidates stay in execution/checkpoints until full147-fit release acceptance.
Do not delete canonical ignored reconstruction payloads during worktree cleanup.
