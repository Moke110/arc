# 06.7 output — tracked record vs local run trees

`output/` is gitignored compute output (the repository's root `.gitignore` ignores
every `output/` tree).  The small evidence files below are force-added with
`git add -f`, exactly like the tracked records of `06.5_tahoe_pretrain/` and
`06.6_cm_random/`.

| path | tracked | content |
|---|---|---|
| `output/init/tahoe_init.json` | yes (`git add -f`) | the three-seed warm-start verification (checkpoint sha256, config digest, epoch, gene axis, loaded key list, empty mismatch lists, anti-random-init canary) |
| `profile/profile_record.json` | yes (`git add -f`) | authoritative bounded wrapper profile: staged checkpoint source/staged sha256, raw `check_init` axis-check output, full warm-start provenance, real 2-epoch run manifest, host/GPU peaks, and the ticket-#6 training-path comparison |
| `units.json` | yes (`git add -f`) | the 45-unit matrix definition + the verified three-seed warm-start block (`sha256`, config digest, gene axis, loaded keys, anti-random-init canary) |
| `record/*` (9 files) | yes (`git add -f`) | candidate record/curves/LR scores/LR selection/15 selected checkpoints, `record.json`, `matrix_run.json`, `resume_check.json`, `record_check.json` |
| `runs/**` | no | the 45 **formal** 06.4 run trees: `manifest.json`, `config.json`, `fit_history.json`, `summary_table.csv`, `run.log`, `checkpoint/model.pt` (weights are gitignored `*.pt`); this tree holds exactly the matrix units |
| `profile/staged_checkpoints/**` | no | the byte-identical snapshot of `#5`'s seed-20261021 checkpoint used by the profile (its sha256 is recorded in the profile record) |
| `profile/runs/**` | no | the bounded profile's one-unit 06.4 run tree (2-epoch checkpoint, manifest, logs); its raw startup lines are embedded in `profile_record.json`, and it is kept out of the formal `runs/` tree |
| `profile/init_seed*.json`, `profile/init_seed*.log`, `profile/profile_queue.json`, `profile/logs/**` | no | the raw axis-check payload/log, the one-unit profile queue and the queue/training logs (their contents are embedded in `profile_record.json`) |
| `input_cache/`, `hash_cache.json`, `queue_status.csv`, `logs/`, `_resume_probe/`, `first_attempt/` | no | derived CSV table cache, frozen-artifact hash memo, queue status, probe tree and the first matrix attempt's evidence (the recorded run is the clean 45/45 re-run after the hash-cache race fix) |

The formal record is **not** written by this worktree: the Lead runs the 45-unit
matrix after ticket #5's sweep finishes, then `make record && make check` writes
and gates `output/record/`.  A profile record is pipeline/throughput evidence
only — it is produced with a 2-epoch bounded budget, never with the locked
100-epoch budget, and `write_record.py --check` refuses a tree without the
completed 45-unit `matrix_run.json`.
