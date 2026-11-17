# Dose-input round1 is an optimization exploration

Status: accepted, 2026-11-17. Amends ADR0075's final-selection, workload,
publication and report-summary obligations following the user's explicit approval
in the round1 grill-with-docs session. Other scientific boundaries remain in force.

Stage 2 03 explores model optimization and does not supply downstream weights.
Round1 comprises the existing 63 cold-drug and 63 held-dose formal fits. Remove
final/deployment selection from active execution and acceptance, and remove model
publication and downstream weight replay from this step's active duties. Keep
cold-drug validation-first method selection and its selected-procedure outer-test
evaluation. Preserve the completed 21 final fits, old model release and original
147-fit acceptance evidence as historical bytes, with an explicit role mapping;
they are excluded from round1 and future dependencies. This avoids making a
method exploration responsible for choosing and maintaining a deployment model.

Report only round1_metrics.png (2×4), round1_learning_curves.png (2×8), and
round1_metrics_cold_drug_dose.png (1×4). Retain plotting data/scripts. Overall
metrics show one test-metric dot per fold×seed run, nine dots per method/scenario;
bars are their equal-weight mean and error bars are mean ± sample SD (ddof=1).
Dot and error-bar colors match bar edges. These are descriptive run spreads,
not confidence intervals or biological replication; held-dose folds also differ
in omitted dose. This explicitly changes plotted bar aggregation from the old
parent-equal OOF summary. Preserve the original scientific OOF summaries and
parent-paired intervals as audit evidence with their distinct interpretation.

Cold-drug dose panels show PCC, Spearman, RMSE and R² versus 0.2/1/10 µM on a
log dose axis, with seven fixed-method curves and the fold-selected procedure.
Use saved parent-equal, equal-seed OOF dose summaries as in the former primary_dose
plot. The learning sheet has cold-drug/held-dose rows, seven native dual-axis
method panels and a method-colored validation-MSE overview in column eight.
Preserve every epoch/run; overview MSE stays logarithmic.

The round1 conclusion is limited to four tested 100-dimensional concatenation
configurations showing no stable prediction benefit. Scaling and FiLM have
different dose mechanisms; no universal architectural claim or test-based
reselection follows. Implementation must reuse existing predictions and metrics,
retain frozen driver/code identities, update active routing and documentation,
and verify a CSV-only reconstruction before replacing report PNGs.

Published specification: [#240](https://github.com/Moke110/arc/issues/240).
Implementation slices: [#241](https://github.com/Moke110/arc/issues/241)
(active scope), then [#242](https://github.com/Moke110/arc/issues/242)
(report acceptance), blocked by #241. Implementation starts in a later session.
