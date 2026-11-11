# log/ — bases execution and audit records

Records for `stage2/01_bases`, per the ADR 0062 step contract: acceptance results,
hashes, resource measurements, migration mappings and compact failure evidence.
Code, tests and disposable intermediates stay in `execution/`; reusable accepted
products stay in `output/`; report figures and their plotted values go to `figure/`.

This directory was created by ticket #155 (parent #154) as a planning and audit
record. No product was moved, rewritten or deleted by #155; every file here
describes the tree as it stood at `dev` `a53c9ff`.

## Records

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

## Reproducing the records

The audit reads metadata and non-expression products only. It never opens
`*.h5ad`, `*.npy`, `*.npz` or `*.bin`, never uses `backed="r"`, and never touches
an external `source/` tree. Measured on this workstation: 0.42 s wall, 130 MiB
peak RSS, 0 expression payloads opened.

```bash
python3 log/155-reproduce.py --step-root /path/to/stage2/01_bases --write /tmp/audit.json
```

`155-reproduce.py` is the audit script text embedded in [`155-audit.md`](155-audit.md);
#156 should promote it to `execution/` (audit code belongs there) instead of
keeping it as record text.

## Conventions

- `*.log` is gitignored repository-wide; records here use `.json`, `.md` and `.csv`.
- Historical sealed records keep their original bytes and their original path
  strings as evidence. Where #155 proposes a path move, the old→new binding is in
  `155-migration-map.json` rather than written into the historical file.
- A record marked `EXPECTED_GAP` in the audit is a documented historical
  dependency on bytes that no longer exist; it is not a current product defect.
