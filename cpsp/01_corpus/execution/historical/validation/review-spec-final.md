## Spec

Targeted rechecks of4ca1cc1 and892abc7 resolve the vehicle-matching and taxonomy findings.

All source adapters now preserve unknown/unconfirmed control matching without treating it as established: every one of6,096 comparisons has `control_matched=False`. Duan additionally records `matched_vehicle_unconfirmed`, retaining the A/R model baseline as source evidence rather than a valid FC denominator. Lightweight current derivation and persisted finalv2 coverage agree: zero computable FCs, zero admitted rows,6,096 exclusions and nine recognized partitions. Scientific exclusions are allowed by#29 and do not establish training readiness.

`GLOSSARY.md` describes five categories and32 active items, distinguishing historical27-item inventory. Release inputs freeze full author-log snapshots, author-log/source-manifest SHA256, adapter hashes and extraction software versions for all five sources.

The earlier three core evidence-boundary findings remain resolved. Immutablev1 remains historical; correctedv2 is the final release under review. Final repository-integrity execution and handoff evidence remain the integrator's responsibility; this targeted check performed no heavy jobs.

Summary:0 outstanding Spec findings.
