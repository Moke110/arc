# CMonge post-training check repair

The first formal CMonge fit completed 10,000 native iterations, then failed a
post-training implementation check. Its selected checkpoint's output was exactly
insensitive to a log-dose +1 probe on the first training condition on CUDA.
Reload error was zero, the Morgan adapter changed, baseline input affected output,
and three of twelve other training conditions responded to dose. The failure
therefore did not establish missing dose wiring or invalid predictions.

The repair makes these single-condition output sensitivities diagnostics for
formal fits. Qualification still requires input participation; finite diagnostics,
changed adapter parameters, exact reload, ordered full-axis finite predictions
and checkpoint identities remain required. No architecture, objective, training
configuration, random draw, dose path or checkpoint-selection rule changed.
Neither test truth nor test scores informed the correction.

## Evidence

- `diagnosis/` retains the checkpoint-only CPU/CUDA reproduction, other-condition
  probes and native artifact-CLI red/green regressions. Harness text snapshots
  preserve the original diagnostic commands; their original output paths are
  historical. The old checkpoint hashes identify the retained binary archives.
- `old_attempts/` retains the original adapter/config/queue/expected identities,
  both attempted runs' compact metadata and logs, guardian telemetry and local
  checkpoint hashes. Model binaries remain local-only.
- `migrate_cmonge_repair.py` is the guarded one-time state migration; it verifies
  the actual host, absence of active jobs, unchanged queue/config/environment
  declarations and 54 non-CMonge identities. Exactly nine CMonge config identities
  change because the adapter source hash changes.
- `migration_dryrun.log` and `migration.json` record checks and execution. The
  migration verifies complete binary archive hashes before clearing old run
  directories for corrected fits. No partial checkpoint is reused as completion.

The user explicitly authorized implementation-bug repair and continuation on
2026-11-06. The failed fit and the second fit stopped gracefully for this repair
remain attempts, rather than disappearing from accounting. Training time carried
forward is 6,810.677218487 seconds for fold0/seed17 and a conservative 1,467.71
seconds for fold0/seed42. Retry counters and batch scale are unchanged. Each run
still has the same 14,400-second cumulative training allowance, with no guarantee
that the corrected fit finishes before its remaining budget.

The saved first checkpoint cannot by itself supply the complete selected-fit
history/checkpoint metadata that the original adapter never wrote after its
exception. Corrected fits therefore restart from the same random seeds with their
prior time charged; no metadata or successful outcome is invented. Both original
binary run directories are retained under the ignored
`output/formal/implementation_repair_archive/cmonge_behavior_gate_v1/` tree.

All 45 completed non-CMonge runs remain reusable. The original serial queue
continues the remaining 18 runs after migration. Scoring remains withheld until
all 63 outcomes settle, including genuine failures/timeouts. The final independent
Standards/Spec review remains pending; the narrow static repair audit found no
P1/P2 issue.
