# tsg_view — do within-condition replicates explain the tenv–Tahoe gap? (first draft; ticket 08 finalizes)

Study step for [ADR 0082](../../docs/adr/0082-view-variants-tenv-tahoe-replicate-structure.md).
It builds expression-view variants from the Bases integer corpus under the ADR 0077
axis-local recipe, freezes each into the Step 04 input format, fits native TranSiGen under
the `no-dose` method and the sealed cold-drug-dose folds, and reports how much of the
tenv–Tahoe Δ gap within-condition replicate structure explains. It publishes comparison
evidence and **no model weights**, and it changes no product of Steps 01–06.

## Views

All views use the frozen G_5000 axis and the ADR 0078 cold-drug-dose folds (split seed 226;
tenv doses 0.2/1/10 µM, `tahoe_bulk` doses 0.05/0.5/5 µM). A condition is
`(cell_line, drug_parent_inchikey, dose_um)`.

| view | rows (tenv / tahoe_bulk) | recipe |
|---|---|---|
| `full` | 4,478 / 66,718 | ADR 0077 unchanged; published as a digest-equal **alias** of the accepted 01 release, rebuilt in staging from corpus integers and asserted against it |
| `merge` | 1,377 / 52,777 | one replicate-pooled profile per condition: sum retained G_5000 integer counts, normalize the sum on its own selected-axis total; control pooled the same way from the union of member control groups |
| `solo` | 1,377 / 52,777 | one treated profile per condition, drawn once per source from a seed frozen before any fit; control is the full-view control |
| `solo-deepest` | 1,377 / 52,777 | as `solo`, choosing the largest retained-count member (sensitivity arm) |
| `full-pooled-label` | 4,478 / 66,718 | full inputs and controls; only each treated profile's Δ becomes its condition's pooled Δ |

The study publishes eight view products (five view families over two sources), freezes them,
and fits seven main views × three folds × training seeds 42/43 plus three `solo-deepest`
fits at seed 42 (45 runs). Evaluation truth is frozen across views.

## Status

Ticket 01 delivers the step skeleton, the shared view harness, the digest-equal `full`
view, the view-aware freeze and verifier, one bounded smoke fit and its tests. The derived
views (tickets 02–04), truth and diagnostics (05), the 45-fit matrix (06), the report and
figures (07) and acceptance (08) are follow-on work; this README is finalized by ticket 08.

## Layout

```
output/views/<view>/<source>.json          # manifests (tracked); payloads ignored
output/runs/<source>/<view>/fold-<f>/seed-<s>/   # fitted runs (ignored)
execution/build.py                         # builder CLI, auto-discovers views/<view>.py
execution/views/__init__.py                # shared harness: readers, normalization, product writer
execution/views/full.py                    # digest-equal full view alias
execution/freeze.py, execution/inputs.py   # view-aware freeze writer and verifier
execution/learner.py, execution/smoke.py   # one-run worker and bounded smoke driver
figure/data/, figure/script/               # plotted values and CSV-only renderer (ticket 07)
log/                                       # build, freeze and acceptance evidence
```

## Reproduction (ticket 01 scope)

Use the `arc` interpreter; the base `python3` cannot reproduce sealed digests. Run one heavy
job at a time with the repository caps.

```bash
PY=/home/chang/miniconda3/envs/arc/bin/python
cd stage2/tsg_view
$PY execution/build.py --view full --source tenv --output output/views \
    --staging execution/checkpoints --repo-root /home/chang/arc --record log/build-01.json
$PY execution/freeze.py --view full --source tenv --views output/views \
    --destination execution/freezes/full/tenv --repo-root /home/chang/arc
$PY execution/smoke.py --view full --source tenv --freeze execution/freezes/full/tenv \
    --output execution/checkpoints/smoke --rows 64 --epochs 2 --device cuda
```

The `full` product carries only the accepted release's manifest digest and payload digests;
the authoritative expression payload remains `stage2/01_bases/output/view-5k/`.
