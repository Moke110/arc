# log/ — bases execution and audit records

Records for `stage2/01_bases`, per the ADR 0062 step contract: acceptance results,
hashes, resource measurements, migration mappings and compact failure evidence.
Code and tests stay in `execution/`; disposable intermediates go to
`execution/checkpoints/`; reusable accepted products stay in `output/`; report
figures and their plotted values go to `figure/`.

## Subdirectories

| path | contents |
|---|---|
| `records/` | former `execution/build_record/`: acceptance, build, replay, launch, resource, failure and publication records (bytes unchanged, including their original absolute path strings). The only addition is the new #156 successor re-issue outputs. |
| `health/final_health/` | former `execution/build_record/final_health/**: repo-health scope, logs, timings and seal ledger. |
| `diagnostic/` | dissolved `output/hvg-metrics-diagnostic-v1/`: the accepted plotted CSV/PNG evidence copies, panel JSON, replay artifacts, per-fold group tables, cohorts evidence, `acceptance.json` and `visual_qa.json`. No expression payload. |
| `historical/shared-v1-mappings/` | the surviving historical v1 mapping CSVs that the sealed identity product binds. |
| `historical/seals/truth-successor-envelope-v1/` | the superseded truth-successor envelope and its sealed R worker, bytes intact, preserved when #156 re-issued the envelope at its stable path. |

## #155 records

| file | contents |
|---|---|
| [`155-inventory.json`](155-inventory.json) | complete machine inventory of every product directory, payload file, seal and execution/acceptance record: size, tracked flag, role, disposition, stable path |
| [`155-inventory.md`](155-inventory.md) | human summary of that inventory: product table, consumer map, acceptance-record index, reconstruction dependencies |
| [`155-audit.md`](155-audit.md) | audit report: what was checked against actually bound metadata, what was reused from prior records, what was newly verified, and the findings that change the plan |
| [`155-metadata-binding-audit.json`](155-metadata-binding-audit.json) | raw audit output: 33 binding checks and 13 structural checks with per-file digests |
| [`155-migration-map.json`](155-migration-map.json) | old path → one stable semantic path, owner, reuse-vs-rebuild, prerequisites, consumer rebinding, acceptance gate and phase |
| [`155-gates.md`](155-gates.md) | reviewable acceptance gates for #156 and #157, the exact tenv truth-fit substitution prerequisites, rejected Tahoe exclusions and the bounded expression checks deferred to #156 |
| [`155-evidence-ledger.json`](155-evidence-ledger.json) | reused-versus-new evidence ledger: per claim, the record that supports it, its digest, and whether #155 re-ran it |
| [`155-health-check.md`](155-health-check.md) | exact scope and outcome of every health check #155 ran, including `check_paths` not completing and why |
| [`155-reproduce.py`](155-reproduce.py) | the #155 audit script text as it ran; its documented successor is `execution/audit_metadata_bindings.py` |

## #156 records

| file | contents |
|---|---|
| [`156-decisions.md`](156-decisions.md) | coordinator decisions D1–D3, the scoped interpretations and their rationale |
| [`156-migration-record.json`](156-migration-record.json) | every old→new rename with per-file inode/size/mtime/sha256 evidence, cleanup, rollback steps and the run command |
| [`156-payload-digests.json`](156-payload-digests.json) | gate E1: one serial re-hash of the accepted payloads (whole-file digests, bytes read, wall time, peak RSS) |
| [`156-successor-reissue.json`](156-successor-reissue.json) | gate E2: candidate create+verify results, superseded and new envelope digests, preserved-bytes path |
| [`156-metadata-binding-audit.json`](156-metadata-binding-audit.json) | the metadata/structure audit re-run against the stable layout |
| [`156-figure-reproduction.json`](156-figure-reproduction.json) | gate E4/E6: complete plotted-values schema, CSV digest, rendered PNG digest, byte-identity result |
| [`156-gates.md`](156-gates.md) | #156 gate outcomes (E1–E5, publication gates, health checks) with commands, environments, resource peaks and the reused-vs-new ledger |
| [`156-evidence-ledger.json`](156-evidence-ledger.json) | machine-readable per-claim ledger: supporting record, its sha256, whether #156 re-ran it, and the outcome |
| [`records/156_truth_successor_{build,verify}.json`](records/156_truth_successor_build.json) | create and independent verify outputs for the re-issued envelope |
| [`records/156_rejected_successor_candidate.json`](records/156_rejected_successor_candidate.json) | the first candidate withdrawn before acceptance (wrong verification-record binding), with its bytes preserved outside `output/` |
| [`records/156_e3_diagnostic_replay.json`](records/156_e3_diagnostic_replay.json) | HVG/diagnostic replay against the stable paths (E3) |
| [`records/156_e5_evaluator_report.json`](records/156_e5_evaluator_report.json) | standalone evaluator report through the re-issued successor (E5) |
| [`health/`](health/) | `156_check_paths.log`, `156_check_tracked.log`, `156_stage2_path_audit.log`, `156_health_script_tests.log` |

## #157 records

| file | contents |
|---|---|
| [`157-final-summary.md`](157-final-summary.md) | the final scientific summary: frozen axis, corpora, 1377-contrast truth with 4206 treated wells and 16596 control memberships, successor split, HVG partitions, the exact fold macro AP values, limitations, provenance, migration pointer, verification scope and the removal record |
| [`157-gates.md`](157-gates.md) | acceptance gates A1–A8 (including the regenerated audit), health checks C1–C3 with exact scope/duration/peak RSS, the final test record, findings F1–F5 and deviations |
| [`157-removals.json`](157-removals.json) | every removed payload/state directory: old path, size, sha256 or tree digest, the rejection/supersession evidence with its own digest, the consumer scan, the reconstruction path, the post-removal check and what was deliberately kept |
| [`157-evidence-ledger.json`](157-evidence-ledger.json) | per-claim reused-versus-new ledger with record digests and the method actually used |
| [`health/157_check_paths.log`](health/157_check_paths.log) | the single allowed `check_paths.sh` run on the branch tree: PASS, 0 unresolved |
| [`health/157_check_tracked.log`](health/157_check_tracked.log) | `check_tracked.sh`: PASS, 242 on disk / 284 tracked / 0 holes |
| [`health/157_verify_all.log`](health/157_verify_all.log) | the single `verify_all.sh` run: 34/34 external items PASS, 4:48.32, 79.4 MiB peak, with the scope distinction |
| [`health/157_tests.log`](health/157_tests.log) | the 13 individual test files with per-file wall time and peak RSS |
| [`health/157_figure_reproduction.log`](health/157_figure_reproduction.log) | CSV-only re-render, byte-identical `35d14266…` |
| [`156-metadata-binding-audit.json`](156-metadata-binding-audit.json) | **regenerated by #157** from the single main checkout: 44 binding + 21 structural checks, `non_pass: []`, one documented `EXPECTED_GAP`, no worktree paths, two consecutive runs byte-identical |

## Reproducing the records

The metadata audit reads manifests, non-expression products, moved records and the
published figure. It never opens `*.h5ad`, `*.npy`, `*.npz` or `*.bin`:

```bash
python3 stage2/01_bases/execution/audit_metadata_bindings.py \
    --step-root /home/chang/arc/stage2/01_bases \
    --write /tmp/audit-after.json
```

`--step-root` is the only override needed for the canonical single-checkout run:
`--figure` and the code root default to that step root, so the record carries no
worktree paths. The audit also reads `log/157-removals.json` (the documented
`EXPECTED_GAP`-style exception for the removed rejected Tahoe payload); with a
split checkout it falls back to the copy next to the script.

The one-off payload re-hash (≈20.9 GB) is
`execution/verify_payload_digests.py`; run it serially through
`execution/capped.sh` with `/usr/bin/time -v`. The figure reproduction is
`execution/publish_diagnostic_figure.py` plus
`figure/script/render_diagnostic_tenv_pr.py`.

## Conventions

- `*.log` is gitignored repository-wide; records here use `.json`, `.md` and `.csv`.
- Historical sealed records keep their original bytes and their original path
  strings as evidence. Old→new bindings live in
  [`156-migration-record.json`](156-migration-record.json) rather than being
  written into the historical files.
- A record marked `EXPECTED_GAP` in the audit is a documented historical
  dependency on bytes that no longer exist; it is not a current product defect.
