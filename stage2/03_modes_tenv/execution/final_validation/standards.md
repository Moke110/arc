# Standards review

Pinned comparison: `git diff f022b1b...cd6d6873085158954b4cd988247acc0fd6c961a6 -- stage2/03_modes_tenv docs/adr/0042-tenv-model-benchmark-redesign.md`. Checked root/Stage 2/step AGENTS, ADR 0030 and 0042, execution-environment rules, and code-review baseline. Parallel Spec review is separate. The supplemental runtime and production overlay/delivery integration received separate-agent source/evidence review. Runtime rereview, compact replay, and final delivery proofs are linked in the coordinator handoff.

## Documented deviations

- **P2, historic joint-v2 resource response.** ADR 0042 lines 100–108 requires five-second observation, withholding dispatch under the RAM floor and stopping the owned task under worsening pressure. `joint_v2_runtime.py:146–158` admits before verifying potentially expensive inputs, and `:192–203` only records active-loop snapshots. The dispatch gap and missing pressure reaction exist in the historical runtime. The actual monitored run had resource headroom. Preserve its evidence and record the implementation limit rather than claim the original runtime enforced live pressure. The supplemental run has separate enforcement.
- **P3, supplemental repo-root discovery.** ADR 0030 lines 67–70 explicitly requires walking to the `.git` directory. `cmonge_supplement.py:30` sets `ROOT = EXECUTION.parents[2]`. It works in the reviewed fixed layout. Retain the exact runtime bytes required by the historical policy digest and record the portability limit. A future successor can use the root-walk idiom already in `cmonge_v2_gpu_adapter.py:8`.

## Judgement calls

- **Duplicated Code.** `joint_v2_apply.py:25–88` and `perturbdiff_v2_apply.py:25–88` repeat the durable publish, inventory and atomic-exchange primitives; `joint_v2_budget.py:23–113` and `perturbdiff_v2_budget.py:23–113` similarly repeat account and process-identity helpers. Divergence made sense across separately reviewed generations, but a future successor could share those small persistence primitives while preserving distinct reviewed contracts.

No violation for the user-specified output figure paths. ADR 0042 lines 145–153 explicitly preserves those figure paths across interim and final publication. Frozen original account files, sources, identities and binary outputs remain separate from the one-slot selection.
