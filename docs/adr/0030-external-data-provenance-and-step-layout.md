# External data provenance and the Stage 2 step layout

**Status: accepted (2026-10-30)**

Two reproducibility gaps have accumulated. First, the four external-source
categories (`pert_data/`, `pert_model/`, `disease_data/`, `drug_emb_model/`) mix
downloaded bytes with locally generated artifacts, and their provenance lives in
five scattered places (`download/manual.md`, `pert_model/STAGE2_DATA.md`, twelve
`VENDORED.txt` files, `pert_data/tahoe/Tahoe_summary.md`,
`pert_data/annotation/README.md`) of which several are not in git at all; 291
symbolic links point at a 40 GB weight archive outside the repository. Second,
`stage2/` steps do not share a layout: code, frozen inputs, run cards, and
per-unit outputs sit in different sibling directories per step, figures sometimes
live two levels down, and the two benchmark substeps each duplicate a full
harness under names (`temp_env_all_18k/5k`) that describe the dataset rather than
the work.

**Decision.**

**1. External item layout (contract one).** Each external category contains one
directory per download source (项), named after the upstream product
(`Tahoe-100M`, `PRnet`, `GSE262419`, `LINCS-LDS-1593`). Each 项 contains exactly
five members:

```
<category>/                      # pert_data | pert_model | disease_data | drug_emb_model
├── AGENTS.md                    # ≤ 30 lines, one line per 项
├── README.md                    # human overview of the category
└── <Item>/                      # one 项 = one download source
    ├── script/                  # downloader + offline verifier (tracked)
    ├── source/                  # downloaded bytes, unmodified, upstream layout (gitignored)
    ├── local/                   # locally generated / locally extended content, declared (gitignored large parts)
    ├── author_log.json          # provenance record (tracked)
    └── manifest.sha256          # per-file hashes over source/ (tracked)
```

`source/` is the untouched download: byte-identical to what upstream shipped, in
upstream's own layout, replaceable only by re-downloading. `local/` is every
addition made here after the download — generated artifacts and locally extended
material — declared one path at a time in `author_log.json`'s
`local_additions[]`. `author_log.json` records the source URL, version, retrieval
time, per-source file count / bytes / tree digest, and the verification result.
`manifest.sha256` lists `source/` only, so `sha256sum -c` verifies a download
offline without touching local material. Downloaders live with the data they
fetch: the Tahoe import scripts move from `stage2/04_tahoe_pretrain/script/` to
`pert_data/Tahoe-100M/script/`, and `fetch_annotation.py` moves to
`pert_data/Annotation/script/`.

**2. Singular `script/`.** Every `scripts/` directory in this repository is named
`script/` from now on.

**3. Stage 2 step layout (contract two).** Every step directory carries
`AGENTS.md` (≤ 30 lines, routing for agents), `README.md` (the one human-facing
summary for that step), `execution/` (harness, configs, run cards, tests, logs,
checkpoints, and other intermediates; agents own this), and `output/` (artifacts
a later step reads or a human must know exists — frozen inputs, predictions,
`.h5ad`). Human-facing `.png`/`.csv` results sit directly in the step directory,
not inside `output/` or `execution/`. A step that splits into
independently-concluded substeps uses third-level `NN.N_<outcome>/` directories
that each carry the full set of those four members plus their own human-facing
figure/table; the second level then holds only those directories plus the step's
`AGENTS.md`/`README.md`. Non-step directories (`_envs/`, `audit/`, `_archive/`)
are exempt and carry a one-paragraph README.

**4. Paths do not encode depth.** Code resolves the repository root by walking up
to the directory containing `.git`, not by counting `parents[N]` or `../..`
segments; the benchmark-root constants are updated once. Frozen evidence files
keep their historical path strings byte-for-byte (a rename must not rewrite
evidence), and a written whitelist governs the old-path audit.

**5. Nothing lives outside the repository.** The 285 weight/data symlinks into
`/mnt/f/arc-weights` are materialized in place (+40 GB); the archive remains as
an offline backup only.

**Confirmed at acceptance (2026-10-30).** The three points spec §10 left open are
settled as follows, and this ADR adopts them as law:

- `local/` **is** allowed as the fifth member of an external 项. It is not an
  escape hatch for sloppiness: it is where declared local additions go, so that
  `source/` stays a pure download that can be diffed directly against upstream.
- The 项 names in spec §2–§5 are adopted verbatim, including their case:
  `Tahoe-100M`, `LINCS-LDS-1593`, `DToxS-iPSC`, `GSE244740`, `GSE247541`,
  `GSE262419`, `XPert-L1000`, `Annotation`; `PRnet`, `CPA`, `chemCPA`, `biolord`,
  `PerturbNet`, `scVIDR`, `Squidiff`, `CMonge`, `DECANT`, `WAVE`, `ExPO`, `State`,
  `TranSiGen`, `XPert`; `GSE108940`, `GSE137293`, `GSE160516`; `HG`, `KPGT`.
- The three human artefacts land as spec §10.3 says: `cm_drugs.html` in
  `pert_data/Annotation/local/`, `SURVEY_2026-10-19_v2.md` in `stage2/audit/`, and
  `01_setup/model_architecture_weights_benchmarks.html` at the step's second
  level (the one non-`.png`/`.csv` exception to the step-level rule).

**Consequences.** Roughly 1,000 path references are rewritten (361
`pert_model/<dir>/`, 73 `pert_data/<dir>/`, 515 `temp_env_all_18k|5k`, 121
depth-relative roots, 11 `04_tahoe_pretrain`), all anchored on full path prefixes
so that dataset identifiers such as `temp_env_all-smiles` and the historical
benchmark id `temp_env_all_5k` survive. `.gitignore` becomes provenance-aware:
`source/` and `local/` stay out of git, while `script/`, `author_log.json`,
`manifest.sha256`, and the item READMEs are tracked; a `check_tracked.sh` guard
asserts that no executable file is silently ignored. Key result tables
(`leaderboard.csv`, `per_drug.csv`, `fidelity.csv`, `run_cards.csv`,
`bootstrap.npz`, `comparison_18k.csv`) are force-tracked, since they are the
conclusions of a 60-unit benchmark. 424 GB of item content is hashed once
(background job) to produce the manifests; the human-opened figures move to
step level; `stage2/PLAN.md` is archived and its live parts absorbed into
`stage2/README.md`. ADR 0029's import paths change; its text stands and gains a
one-line note. Materializing 40 GB of weights makes the repository self-contained
at the cost of disk, accepted because the host has 289 GB free.

See `stage2/README.md`, `pert_data/AGENTS.md`, `pert_model/AGENTS.md`,
`.scratch/data-provenance-restructure/spec.md`. These are final destinations:
`stage2/README.md` does not exist until ticket 08, and `pert_data/AGENTS.md` /
`pert_model/AGENTS.md` reach their final one-line-per-项 form in tickets 03 and 04.

## Clarification recorded at acceptance (2026-10-30)

Two readings needed spelling out once the tree was built; neither changes a decision.

1. **§2 `script/` singular applies to our own directories.** A 项's `source/` tree is
   byte-identical upstream content and keeps upstream's own layout, so
   `drug_emb_model/KPGT/source/scripts/` stays plural — renaming it would break the
   "diff against a fresh download" property `source/` exists for. `stage1/` is likewise
   out of scope (§1 rule 5) and its frozen `xpert-fda-miri-v0/scripts/` is untouched.
2. **§4's "paths do not encode depth" is a Stage 2 rule.** The finalized Stage 1
   workflows keep their existing `Path.parents` arithmetic; they were not restructured.

Everything else the ADR states is what the tree implements: 27 项 across the four
categories, `AGENTS.md`/`README.md` at every category and every Stage 2 step (each
`AGENTS.md` <= 30 lines), four Stage 2 steps with `03_models_cmbench` split into
`03.1_cold_drug_18k` and `03.2_hvg_sensitivity_5k`, every `source/` verifiable offline
with `sha256sum -c`, and no content outside the repository.
