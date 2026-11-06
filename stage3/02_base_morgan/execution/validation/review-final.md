# Final two-axis review — #64

Reviewed feature changes from `c17686dcf1409d8b069c8540685ddf28a3a9a1e3`
through `2fd12768b755afe405f14a41f726e484c2f9812d`:

```bash
git diff c17686d...2fd1276 -- AGENTS.md stage3/AGENTS.md stage3/README.md stage3/02_base_morgan
```

Unrelated interleaved Stage 2 work is outside this review. Separate fresh-context
reviewers inspected code, tests, routing and release evidence. The large split
manifest was inspected through compact schema/count/hash queries, not pasted into
review context. Neither reviewer ran tests, training or heavyweight imports while
the repository provenance sweep was active.

## Standards

No actionable findings. Root/Stage 3 routing, ADRs 0030/0045 and the execution
environment rules are satisfied. Sequential dispatch, explicit thread limits,
resource admission, provenance, immutable completed results and preserved failed
attempts follow repository conventions. The documented portable memory-limit
fallback permits RLIMIT_AS. Small repetitions did not warrant refactoring the
frozen protocol.

## Spec

No actionable findings. Features, architecture, optimizer, seeds, validation-only
selection, training-only dose scaling, three partition settings, metrics,
provenance and checkpoint replay match #58. There are 60 unique completed formal
runs: 18 cold-drug, 24 cold-dose and 18 warm-start, with 300 persisted run files
and no failed/pending runs. The final summary hash matches the independent
numerical audit. Full completion requires the canonical 60-fit matrix. Source
normalization and biological/calibration limitations are explicit.

Pre-review findings were resolved in `ae2426a`: single-thread sparse cold-dose
solving and separate selected-scope/full-protocol completion. Resource admission
was refined in `a33ffe1` after a small maintenance job caused a conservative
pause. Training code, data, partitions and model configuration remained frozen.

Final findings: Standards 0; Spec 0. The independently verified formal data record
is [formal-data-audit.md](formal-data-audit.md). Repository-wide provenance
verification is a separate **unfinished external-resource gate**, recorded in
[offline-verification.md](offline-verification.md); it is not claimed to pass.
