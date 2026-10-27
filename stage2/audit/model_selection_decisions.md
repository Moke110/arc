# Stage 2 — current model-selection decisions

This record governs follow-up model selection. Historical benchmark registries,
leaderboards, and entrant audits describe the experiments actually performed;
they do not override the decisions below.

## 2026-10-26 — exclude Squidiff from further selection

**Status: accepted by the user.** The delivery requirement remains prediction
of 15,000+ genes. Squidiff is excluded from subsequent Stage 2 model selection,
tuning, and MIRI adaptation. Do not allocate further training or debugging runs
to it unless the user explicitly reopens this decision.

**Reason:** Squidiff failed to converge under the executed 18,164-gene and
5,000-gene recipes. Continuing to diagnose and adapt it has uncertain benefit
and would divert effort from the remaining candidates. This is a project
resource-allocation decision, not evidence that the method is intrinsically
incapable of whole-transcriptome prediction.

The exclusion is **not** based on a 200-gene architectural limit or an assumed
order-of-magnitude parameter increase. The architecture is configurable; the
18k model has 164,609,840 parameters (about 3.04 times the 200-gene model), and
its executed five-fold run took about 4.71 hours, within the previous budget.
Successful execution did not produce an adequately trained model.

**Remaining candidates:** PRnet, CPA, chemCPA, biolord, PerturbNet, and CMonge.
Their presence here does not establish convergence, fidelity, or final winner
status. Keep the CPU baselines, especially `delta_mean`, in follow-up comparisons.

**Preservation:** retain Squidiff's vendored code, environments, setup records,
configs, checkpoints, predictions, and benchmark rows. Leave both archived
seven-model benchmark registries and their 60-unit results intact for
reproduction. New comparison queues must use the six candidates above.

Evidence:

- [Fairness and resource audit](benchmark_fairness_2026-10-26.md)
- [18k report](../03_bench_models/temp_env_all_18k/REPORT.md)
- [5k report](../03_bench_models/temp_env_all_5k/REPORT.md)
- [18k Squidiff run card](../03_bench_models/temp_env_all_18k/models/squidiff/README.md)
