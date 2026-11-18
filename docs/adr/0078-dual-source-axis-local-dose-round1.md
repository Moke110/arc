# Rebuild dose-input round1 on two axis-local 5k views

**Status** — Accepted design, November 18, 2026. Publication follows full candidate acceptance.

Apply the paired ADR 0077 views to independent Tahoe and tenv comparisons.
Freeze source, profile and matched-control identities and the view recipe,
normalization, axis and payload digests. Use all 65,318 Tahoe and 4,206 tenv
treated samples, preserving their real dose units and cell contexts.

Amend ADR 0075 and ADR 0076 for the successor generation. Retain seven fixed
dose methods with 32-dimensional dose concatenations, three outer folds and
training seed 42 per source and scenario. Run 84 formal fits. Keep the native
backbone and strict best-validation treated-expression MSE checkpoint
restoration. In cold-drug, hold out parents across all doses and cell lines.
In cold-drug-dose, pair each held parent group with one omitted dose and test
only their intersection. Fit and validate using other parents at the other
doses. Record both excluded strips. Use split seed 226 and hold out the ceiling
of 10% of development parents for validation.

Evaluate Delta per treated sample against its own matched control, then pool
sample metrics equally. Keep constant and undefined counts. For cold-drug,
lock method choices from validation Delta PCC before testing and reuse the
locked checkpoints for selected-procedure evaluation. Cold-drug-dose compares
fixed methods. Its test population covers only the held intersections.
Descriptive equal-fold mean ± sample SD is reported separately from pooled
sample statistics.

Publish three CSV-reproducible figures with Tahoe above tenv, visible zero
references and explicit axis breaks where needed. Replace active report
products after exact fresh-process replay, scientific and report acceptance,
and the offline health checks. Retain the historical tenv input, fit, weight,
report and scientific identities. There is no final fit or downstream weight
publication. Stage 2 02 and Stage 3 remain independent.

- [ADR 0075](0075-transigen-all-dose-input-selection.md)
- [ADR 0076](0076-dose-input-round1-exploration-boundary.md)
- [ADR 0077](0077-five-k-only-axis-local-normalization.md)
- [Specification #246](https://github.com/Moke110/arc/issues/246)
