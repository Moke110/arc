# Stage 2 step path migration (#268)

ADR 0079 inserts the drug representation comparison before dose and transfer.
Resolve a sealed path by substituting its exact step-root prefix and keeping the
remaining components. For an absolute path pinned to an earlier checkout, first
interpret the suffix beginning with that checkout's `stage2/` relative to the
current repository. A historical external or temporary location has its original
meaning and requires the recovery instructions for that run.

| Original root in evidence | Current root |
| --- | --- |
| `stage2/03_tsg_dose/` | `stage2/04_tsg_dose/` |
| `stage2/04_tsg_transfer/` | `stage2/05_tsg_transfer/` |

Apply the same mapping to ignored `execution/checkpoints/` and dose
`output/model/`. The [pre-move inventory](pre-move-ignored-inventory.json)
records each payload group in the workstation checkout before relocation.
It excludes Python bytecode caches. The inventory digest is SHA-256 of the
compact JSON encoding of the sorted triples
`[path_relative_to_group, lstat_size, symlink_target_or_null]`. The group root
is `.` when the group is a single file. The digest verifies the inventory
of names, sizes and links. Per-file scientific integrity stays bound to the
original receipt hashes.

Move the ignored subtrees on the same filesystem, retaining every relative
path. Rebuildable caches may be removed separately. Move the subtrees out of
the old roots before switching branches when doing so avoids any worktree
path collision. After checkout, place them at the corresponding new roots
and verify the inventory and the documented downstream reads before running
the Stage 2 path gate. Data payloads are checkout-local.

Sealed receipts and the recorded paths in them are retained byte for byte.
When executing an old plan whose original checkout path is part of its
contract, use a checkout of the original revision with its expected paths
and the recovered inputs. Current readers and report entry points use the
new locations.
