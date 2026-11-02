# 06.5 output — tracked record vs local run trees

`output/` is gitignored compute output (the repository's root `.gitignore` ignores
every `output/` tree).  The small evidence files below are force-added with
`git add -f`, exactly like the tracked records of `06.1_eligibility/`,
`06.4_pipeline/output/smoke/` and the planned `06.8_baselines/` record.

| path | tracked | content |
|---|---|---|
| `profile/profile_record.json` | yes (`git add -f`) | authoritative short-profile record: input identity, resolved budget per unit, validation curve, wall/epochs/steps, host and GPU peaks, stop reason |
| `profile/profile_record_curve.csv` | yes (`git add -f`) | one row per (unit, epoch) of the profile |
| `profile/model_cards/<unit_id>/model_card.md` | yes (`git add -f`) | the 06.4 pipeline model card per profile unit |
| `formal/tahoe_pretrain_record.json` | yes, after the Lead-run sweep | the three-seed formal record (checkpoints' path + sha256, config digests, full curves, resources) |
| `formal/tahoe_pretrain_curve.csv` | yes, after the sweep | the three validation AUPRC@200 curves |
| `formal/model_cards/tahoe_seed*/model_card.md` | yes, after the sweep | one model card per seed |
| `runs/**` | no | the 06.4 run trees: `manifest.json`, `config.json`, `fit_history.json`, `summary_table.csv`, `run.log`, `checkpoint/model.pt` (weights are gitignored `*.pt`) |
| `profile/runs/**`, `formal/runs/**/hash_cache.json` | no | local run trees and derived hash memos |
| `formal/tahoe_pretrain_record_partial.json` | no | a per-seed pre-sweep record; it is **not** the formal evidence and says `complete: false` |
| `*/06.5_unit_resources.json` | no | the wrapper's sampled GPU/RSS numbers, kept beside each run for `--resume` |

The record JSON is the authoritative artifact; the CSV is a convenience view of
its curves.  A checkpoint is identified by path + sha256 inside the record and
stays local (regenerable from the seed and the frozen input).

A **profile** record is pipeline/throughput evidence only: it is produced with
the 2-epoch bounded `--smoke` budget and/or one full-width epoch, never with the
formal 100-epoch locked budget.  `check_record.py --expect-mode formal` refuses
it.
