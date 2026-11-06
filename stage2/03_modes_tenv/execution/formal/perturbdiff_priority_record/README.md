# PerturbDiff dispatch priority

The user authorized an order-only change on 2026-11-06: finish the active
CMonge fold 0 seed 42 naturally, dispatch all nine PerturbDiff runs, then the
seven remaining CMonge runs. `queue.before.json` preserves the original queue;
`order_audit.json` proves the same 63 job objects and all other queue fields.
Frozen identities, commands, model configurations, budgets and adapters are
unchanged. The original queue path remains the continuation entry point.

While this change was being prepared, actual-host inspection found dispatcher
938762 absent and its worker 972642 still running. No controller had been
launched and no signal had been sent by the priority implementer. The cause of
the dispatcher's disappearance is not established. `emergency_guardian.py` was
launched as PID 1000033 to restore five-second resource telemetry, the original
pressure policy, cumulative four-hour training accounting (including the
unsampled gap), and the separate inference watchdog. The scope's memory/swap
limits remained in place. Independent review verified the live guardian.

`drain_boundary.py` waits for the original worker and guardian to finish, then
acquires the scheduler lock, rejects surviving scope members or other running
states, preserves guardian forced outcomes, and validates successful artifacts
in a subprocess with a 300-second limit. Success uses native checkpoint times
plus the existing 1467.71-second training carry. It preserves attempts and retry
accounting before executing the original scheduler against the reordered queue.
It performs no test scoring. Runtime JSON evidence records readiness and the
actual natural boundary. The finalizer/scheduler runs as a durable user service,
`arc-tenv-priority-boundary.service`, independent of tool-session retirement.
