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

## Reproducing the records

The metadata audit reads manifests, non-expression products, moved records and the
published figure. It never opens `*.h5ad`, `*.npy`, `*.npz` or `*.bin`:

```bash
python3 stage2/01_bases/execution/audit_metadata_bindings.py \
    --step-root /home/chang/arc/stage2/01_bases \
    --figure /home/chang/arc/stage2/01_bases/figure \
    --write /tmp/audit-after.json
```

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
