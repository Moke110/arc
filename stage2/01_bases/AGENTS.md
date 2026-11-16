# AGENTS.md — stage2/01_bases

Routing for the accepted Stage 2 foundation. Read this, then the owning section of
[`README.md`](README.md); commands and contracts live there.

## Read first
- [`README.md`](README.md) — per-product contracts, build/replay commands, layout,
  readers, errata and removals.
- [`log/acceptance.md`](log/acceptance.md) — acceptance replays, digests, resource
  peaks, the old→new mapping and "no retired path read" evidence.
- [`log/README.md`](log/README.md) — the `log/*_verify.json` records and removal record.
- ADRs: 0072 (uniform disease gate and unified TOX successor), 0071 (the current gene-axis contract — DEG-only tiered sources,
  `DEG2000`/`f2000`), 0066 (the amended contract), 0068 (fusion amendment kept),
  0067 (errata), 0069/0070 (`AF` removals, superseded by 0071), 0062 (products),
  0061 (frozen axis), 0030 (external-data and step layout).

## Products
`output/raw_counts/` → `output/corpus/` → `output/view-2k/`, `output/view-5k/`;
`output/disease-deg/**` (five declared successor categories; legacy output keeps its manifest roster until promotion),
`output/disease-fusion/` (the ADR 0071 per-disease level-1 rank + `DEG2000`),
`output/gene-axis/` (`G_2000`/`G_5000`, the `DEG2000` coverage table and the
human-facing `gene_axis_fill_rule.docx`), `output/gene-identity/`, `output/tenv-hvg/`,
`output/fold_n_seed/`, `output/palette/` and `output/api/` (the single live public
implementation).

**One v1 filename is deliberate.** `output/gene-axis/source_topn_coverage.csv`
keeps its v1 *name* but carries v2 *semantics*: long format
(`N, n_axis, group, group_label, numerator, denominator, rate, note`) with each
group's own denominator, i.e. the `DEG2000` coverage table, not a top-`N` table.
Renaming it would re-pin the gene-axis manifest, both view manifests and force a
republish, so the name is retained for compatibility — do not "fix" it.

## Rules
- Bind `output/api/bases/*` through `execution/public_api.py`; never copy a module.
- Drivers are `publish|verify`, staging candidates under `execution/checkpoints/`.
- R interpreters resolve through `output/api/bases/r_env.py` (`ARC_RSCRIPT` /
  `ARC_R_HOME`; see README).
- Readers only; never open a payload h5ad with `backed="r"` or materialize a dense observation × 18,164 matrix.
- Heavy work: one job at a time through `execution/capped.sh` with `/usr/bin/time -v`;
  run `script/verify_all.sh` alone.
- `log/` is execution evidence, not a change history — git is the history.
