# Stage 2 historical replay

The original ADR 0082 study is preserved in [tsg_view](tsg_view/README.md).
Its 45-fit report, original recipes, and release and fit receipts keep their
original identities. The snapshots below reproduce the checkout paths used
by the seals.

| Original path | Current payload location |
| --- | --- |
| `stage2/01_bases` | `stage2/01_bases`, recovery and shared external derivation kernels |
| `stage2/02_models_tenv-5k` | `stage2/04_models_tenv-5k` |
| `stage2/03_tsg_drug_emb` | `stage2/05_tsg_drug_emb` |
| `stage2/03_tsg_drug_emb/output/representations` | `stage2/util/output/representations` |
| `stage2/04_tsg_dose` | `stage2/06_tsg_dose` |
| `stage2/05_tsg_width_dose` | `stage2/07_tsg_width_dose` |
| `stage2/06_tsg_transfer` | `stage2/08_tsg_transfer` |
| `stage2/tsg_view` | `stage2/03_view/03.1_hvg-deg_views/historical/tsg_view` |

## Move the materialized checkout once

Run the helper from the candidate checkout before merging it into the
materialized checkout. Store the plan outside Git. It records ignored Stage 2
files under the mapping by inode, size and modification time. Quiesce Stage 2
writers while making and applying the plan.

```sh
CANDIDATE=/path/to/migrated/checkout
MATERIALIZED=/path/to/materialized/checkout
PLAN=/path/to/stage2-payload-migration.json
export PYTHONDONTWRITEBYTECODE=1
python "$CANDIDATE/script/migrate_stage2_payloads.py" plan \
  --root "$MATERIALIZED" --plan "$PLAN"
```

Before switching the materialized checkout, `link-active` can expose existing
payloads at the new names in the candidate checkout for read-only validation.
Use `--target "$CANDIDATE"` with the same root and plan. Expose required
external sources read-only and verify them against their source manifests.
Write candidate outputs to separate directories.

After the Git merge, run the helper with `apply`. It preflights collisions,
unplanned files and identity drift before renaming anything. It retains file
identities and rolls back the current batch on an error. Keep the plan for
replay or use `reverse` before restoring the old Git layout. Publish new
corpus, view and input artifacts to their versioned paths separately.

## Replay with the original code

Use commit `4080ce9187906ae9b0193e26414f125de45fd09d` for the original five
consumer steps. Use commit `212d3e0` for the ADR 0082 study and its Steps 01–06
inventory, which predates the later gene-axis presentation update. Both
commits use the original path names.

```sh
SNAPSHOT=/path/to/replay/worktree
git -C "$MATERIALIZED" worktree add --detach "$SNAPSHOT" 212d3e0
python "$CANDIDATE/script/migrate_stage2_payloads.py" link-snapshot \
  --root "$MATERIALIZED" --plan "$PLAN" --target "$SNAPSHOT"
```

Replay links retain the original names and physical file identities on the
same filesystem. Follow the snapshot's owning AGENTS.md. Run read-only checks against linked
files and write figures and receipts to separate candidate paths. The ADR 0082
receipt uses the original compact `.scratch/tsg-view/research/baseline`
inventory and planning ref. Supply those with `--baseline` and `--refs`, or
place the original inventory in the snapshot. Verify any external sources the
historical workflow names against their manifests.
