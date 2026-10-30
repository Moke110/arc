# AGENTS.md — stage2/04_tahoe_data

Routing index. Purpose, corpus numbers and reproduce commands live in
`README.md`; the contract is ADR 0029. This step is **post-download processing
only**: it downloads nothing and imports nothing.

```
04_tahoe_data/
├── AGENTS.md  README.md
├── execution/build_pseudobulk.py     scan -> aggregate -> assemble -> qa
└── output/pseudobulk/                the corpus and its QA record
```

- `execution/build_pseudobulk.py` — the only entry point. `--check` prints the
  resolved repository root, the Tahoe input and the output directory. It finds
  the root by walking up to `.git`; never add `parents[N]`/`../..` arithmetic.
- Input is the imported atlas `pert_data/Tahoe-100M/source/`; the importer and
  its verifiers live in `pert_data/Tahoe-100M/script/` (`fetch_tahoe.py`,
  `tahoe-dl.sh`, `run_fetch.sh`, `verify_import.py`, offline `verify.sh`).
- `output/pseudobulk/tahoe_pseudobulk.h5ad` — 25.8 GB corpus, local-only, gitignored.
- `output/pseudobulk/{qa.json,qa_report.md,qa_contrasts.csv,manifest.csv}` —
  tracked build record: corpus shape, QA gates, per-contrast agreement, profile table.
- `output/pseudobulk/_parts/` — 14 per-plate CSR parts, the resumable part cache
  (intermediate, gitignored; `aggregate` rewrites a missing part).
- `output/pseudobulk/_state/` — `shard_plate.json` scan bookkeeping (machine
  state, gitignored; safe to delete, `scan` rebuilds it).

From the repository root:

```bash
python stage2/04_tahoe_data/execution/build_pseudobulk.py --check
python stage2/04_tahoe_data/execution/build_pseudobulk.py qa      # rebuild only the QA
```

`--plates plate4,plate5` limits any stage; `--force` rebuilds existing parts.
