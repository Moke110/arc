# AGENTS.md — stage2/01_bases

Routing for the accepted Stage 2 foundation. Read this, then the owning section of
[`README.md`](README.md); commands and contracts live there.

## Read first
- [`README.md`](README.md) — per-product contracts, build/replay commands, layout,
  readers, errata and removals.
- [`log/acceptance.md`](log/acceptance.md) — acceptance replays, digests, resource
  peaks, the old→new mapping and "no retired path read" evidence.
- [`log/README.md`](log/README.md) — the `log/*_verify.json` records and removal record.
- ADRs: 0066 (contract), 0068 (fusion amendment), 0067 (errata), 0069 (AF removal),
  0062 (products), 0061 (frozen axis), 0030 (external-data and step layout).

## Products
`output/raw_counts/` → `output/corpus/` → `output/view-2k/`, `output/view-5k/`;
`output/disease-deg/**` (five categories + the `cardiotox/` crosswalk),
`output/disease-fusion/` (the ADR 0071 per-disease level-1 rank + `DEG2000`),
`output/gene-axis/` (`G_2000`/`G_5000`, the source top-`N` coverage table and the
human-facing `gene_axis_fill_rule.docx`), `output/gene-identity/`, `output/tenv-hvg/`,
`output/fold_n_seed/`, `output/palette/` and `output/api/` (the single live public
implementation).

## Rules
- Bind `output/api/bases/*` through `execution/public_api.py`; never copy a module.
- Drivers are `publish|verify`, staging candidates under `execution/checkpoints/`.
- R interpreters resolve through `output/api/bases/r_env.py` (`ARC_RSCRIPT` /
  `ARC_R_HOME`; see README).
- Readers only; never open a payload h5ad with `backed="r"` or materialize a dense observation × 18,164 matrix.
- Heavy work: one job at a time through `execution/capped.sh` with `/usr/bin/time -v`;
  run `script/verify_all.sh` alone.
- `log/` is execution evidence, not a change history — git is the history.
