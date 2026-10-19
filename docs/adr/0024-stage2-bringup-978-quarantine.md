# Stage 2 bring-up quarantines the 978-landmark PRnet demo

**Status: accepted (2026-10-18)**

ADR 0020 requires Stage 2 to use no 978-landmark subset anywhere. The PRnet
public checkpoint that ships with the repository is a 978-gene LINCS model, and
its native `test_demo.py` is the fastest way to prove that the PRnet code and
weights work on this machine.

**Decision.** Split Stage 2 into a **bring-up** stage and a later **selection**
stage. The 978-gene PRnet demo is allowed **only** inside
`stage2/01_setup/prnet/`, and its outputs are labelled
`bring-up-only`. No Stage 2 evaluation, metric, fold, or model comparison may
read a 978 artifact.

**Consequences.** Bring-up can use each model's native asset without pretending
it is a Stage 2 result. PRnet's run card records the quarantine. The
whole-transcriptome rule of ADR 0020 still governs the selection stage. The
bring-up outputs are not comparable with the other models and must not be
reported as a leaderboard.

See `stage2/PLAN.md` and `stage2/01_setup/prnet/README.md`.

**Update (2026-10-19, ADR 0025).** `01_bringup/` was renamed `01_setup/`; the
active label is **setup-only**. “Bring-up” is retained above as the historical
term for this decision.
