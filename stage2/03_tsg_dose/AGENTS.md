# TranSiGen dose-input optimization — round1

Read README.md, ../AGENTS.md, ADR0075 and its ADR0076 amendment before changes.
Active planning, training and acceptance cover63 cold-drug plus63 held-dose fits.
Use execution/run_matrix.py and execution/transigen.py; both exclude final fits.
This step supplies comparison evidence and three round1 report PNGs, with no
active downstream weight contract. Model publication is retired.

Public scientific implementation stays in output/api/tsg_dose. Its accepted
bytes and the frozen historical drivers in execution/historical preserve original
fit identities; log/round1-role-mapping.json resolves those source hashes.
Keep the21 completed final fits, old output/model, original147-fit registration
and sealed acceptance records unchanged as historical reconstruction evidence.
Read inputs through bounded Bases readers; source bytes, accepted02 and Stage3
scientific artifacts stay unchanged.

Tests use the approved input-freeze, learner/predict/replay, plan/execute and
score/accept/report boundaries. Preserve strict best-checkpoint restoration,
validation-first cold-drug choice locks and selected-procedure outer testing.
Held-dose never selects a winner. Retained126 origins require exact replay proof;
round1 acceptance does not require final origins or downstream weight replay.

Report changes reuse saved tables/predictions without training or re-prediction.
Use execution/publish_report.py prepare/accept/install for report-only
publication. Build and verify candidates, reproduce the three PNGs from CSVs
alone and inspect them before replacing figure outputs.
The active CSV renderer emits exactly round1_metrics.png,
round1_learning_curves.png and round1_metrics_cold_drug_dose.png.
Historical12-PNG reconstruction uses execution/historical/render_v1.py. Retain original scientific summaries with
their parent-equal interpretation. Formal fits, if separately requested, run singly
under capped.sh and /usr/bin/time -v with memory preflight and10-second RAM/GPU/WSL
monitoring. Preserve canonical ignored reconstruction payloads during cleanup.


ADR0078 dual-source release candidates follow the release workflow in README.md.
Keep the sealed 126-fit tree active until installation. After installation,
route the legacy CSV renderer to `execution/historical/report_126/`. The
input preparer bytes remain bound to the current ADR0078 input freeze.
