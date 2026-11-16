# AGENTS.md — stage2/01_bases

Routing for the current ADR 0073 foundation; read the owning README section.

## Read first
- [`README.md`](README.md): live contracts, public readers, build/replay commands.
- [`log/acceptance.md`](log/acceptance.md): current receipt and sealed predecessors.
- [`log/README.md`](log/README.md): replay, document and source-proof records.
- ADR 0073: all-source collection, tier-first DEG1000/f1000, four nested axes.
  ADRs 0071/0072 retain historical contracts; 0062 defines products, 0061 identity,
  0030 external provenance. Apply `docs/agents/execution-environment.md`.

## Products and rules
- `raw_counts/` → `corpus/` → `output/view-{2,3,4,5}k/` (tenv/tahoe_bulk each).
- `output/disease-deg/`: five categories; `disease-fusion/`: full union + DEG1000.
- `output/gene-axis/`: four axes, candidate/composition/coverage tables + manifest.
- `gene_axis_fill_rule.docx`: step-root human document; figures/CSVs under `figure/`.
- `source_topn_coverage.csv` retains its filename but uses explicit DEG1000/group denominators.
- Bind the single `output/api/bases/*` implementation through `execution/public_api.py`.
- Drivers stage `publish|verify` candidates under `execution/checkpoints/`, then replace.
- Resolve R through `bases.r_env`; README distinguishes registry R from TOX R.
- Read payloads through bounded public readers; never open h5ad with `backed="r"`.
- Run one heavy job at a time with `execution/capped.sh` and `/usr/bin/time -v`;
  run `script/verify_all.sh` alone. Preserve ignored payloads before worktree cleanup.
- `log/` is execution evidence; sealed historical records keep their original bytes.
