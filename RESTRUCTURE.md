# Repository restructure — external-data provenance and the Stage 2 step layout

**Status:** implemented on branch `feat/data-provenance-restructure`, 17 commits ahead of
`main` (`ae61965`), working tree clean, not yet merged.
**Date:** executed 2026-10-30; this summary written 2026-10-31.
**Decision record:** [ADR 0030](docs/adr/0030-external-data-provenance-and-step-layout.md).
**Scope:** 781 files changed, +99,761 / −44,367 against `main`.

The repository had two reproducibility gaps. The four external-source categories mixed
downloaded bytes with locally generated artifacts, and their provenance lived in five
scattered places of which several were not in git at all; 285 symbolic links pointed at a
40 GB weight archive outside the repository. Separately, Stage 2 steps did not share a
layout: code, frozen inputs, run cards and per-unit outputs sat in different sibling
directories per step, figures sometimes lived two levels down, and the two benchmark
substeps each duplicated a harness under names that described the dataset rather than the
work.

This restructure closes both. Every external source is now one **项** (item) with its own
downloader, untouched source tree, declared local additions, provenance record and offline
hash manifest. Every Stage 2 step now has the same four members plus two documents. Nothing
the repository references lives outside it, and every claim is checkable with three scripts
that need no network.

---

## 1. The two contracts

### 1.1 External item (项)

Each external category contains one directory per download source:

```
<category>/                      # pert_data | pert_model | disease_data | drug_emb_model
├── AGENTS.md                    # <= 30 lines, one line per 项
├── README.md                    # human overview of the category
└── <Item>/
    ├── script/                  # downloader + offline verifier            [tracked]
    ├── source/                  # untouched download, upstream layout      [ignored]
    ├── local/                   # declared local additions (only where they exist) [ignored]
    ├── author_log.json          # provenance record                       [tracked]
    └── manifest.sha256          # per-file hashes over source/ only       [tracked]
```

Rules that make it auditable:

- **`source/` is byte-identical to upstream.** Never edited in place, replaceable only by
  re-downloading. That is what lets it be diffed against a fresh download.
- **`local/` is every addition made here**, declared one path at a time in
  `author_log.json` → `local_additions[]`.
- **`manifest.sha256` covers `source/` only**, and `script/verify.sh` is
  `cd source && sha256sum -c ../manifest.sha256`. A regenerating local artefact can
  therefore never fail the download verification.
- **`tree_sha256`** is `sha256` over the manifest file bytes, where the manifest is sorted
  by source-relative path.
- **Every `scripts/` directory in our own trees is named `script/`.** A vendored `source/`
  tree keeps upstream's own layout, so `drug_emb_model/KPGT/source/scripts/` stays plural.

### 1.2 Stage 2 step

```
stage2/<NN>_<name>/
├── AGENTS.md        # <= 30 lines, routing only
├── README.md        # the one human summary for this step
├── execution/       # harness, configs, run cards, tests, logs, intermediates (agents own)
├── output/          # artifacts a later step reads or a human must know exists
└── *.png / *.csv    # human-facing conclusions sit at step level
```

A step with independently concluded parts splits into `NN.N_<outcome>/` substeps that each
carry the full set; the parent then holds only its two documents and the substep
directories. Non-step directories (`_envs/`, `audit/`, `_archive/`) are exempt and carry a
one-paragraph README.

Two binding consequences: **code never finds the repository root by counting parents** — it
walks up until it sees `.git` and asserts it got there; and **no referenced content lives
outside the repository**.

---

## 2. External data: 27 items, 453 GB, verifiable offline

Eight in `pert_data/`, fourteen in `pert_model/`, three in `disease_data/`, two in
`drug_emb_model/`. Sizes are the measured bytes under each `source/` (6,005 files total):

| 项 | files | bytes | size |
|---|---:|---:|---:|
| `pert_data/Tahoe-100M` | 4,430 | 428,811,080,261 | 428.81 GB |
| `pert_data/XPert-L1000` | 3 | 9,833,358,088 | 9.83 GB |
| `pert_data/DToxS-iPSC` | 17 | 2,657,111,110 | 2.66 GB |
| `pert_data/LINCS-LDS-1593` | 213 | 913,276,012 | 0.91 GB |
| `pert_data/GSE262419` | 18 | 83,444,654 | 0.08 GB |
| `pert_data/GSE244740` | 2 | 37,480,421 | 0.04 GB |
| `pert_data/Annotation` | 1 | 16,982,321 | 0.02 GB |
| `pert_data/GSE247541` | 2 | 999,568 | 0.00 GB |
| `pert_model/State` | 162 | 4,356,241,170 | 4.36 GB |
| `pert_model/PerturbNet` | 138 | 2,206,310,646 | 2.21 GB |
| `pert_model/CPA` | 50 | 1,534,453,059 | 1.53 GB |
| `pert_model/chemCPA` | 271 | 1,059,093,253 | 1.06 GB |
| `pert_model/scVIDR` | 273 | 165,168,207 | 0.17 GB |
| `pert_model/XPert` | 14 | 155,247 | 0.00 GB |
| `pert_model/biolord` | 43 | 180,822,754 | 0.18 GB |
| `pert_model/PRnet` | 28 | 109,148,389 | 0.11 GB |
| `pert_model/TranSiGen` | 45 | 67,642,560 | 0.07 GB |
| `pert_model/CMonge` | 85 | 34,029,021 | 0.03 GB |
| `pert_model/WAVE` | 24 | 1,523,907 | 0.00 GB |
| `pert_model/Squidiff` | 26 | 1,234,342 | 0.00 GB |
| `pert_model/DECANT` | 17 | 556,363 | 0.00 GB |
| `pert_model/ExPO` | 83 | 212,744 | 0.00 GB |
| `drug_emb_model/KPGT` | 44 | 929,958,685 | 0.93 GB |
| `disease_data/GSE108940` | 1 | 2,766,695 | 0.00 GB |
| `disease_data/GSE137293` | 2 | 2,582,175 | 0.00 GB |
| `disease_data/GSE160516` | 1 | 2,004,508 | 0.00 GB |
| `drug_emb_model/HG` | 12 | 136,317 | 0.00 GB |


**Total: 6,005 files / 453,007,772,477 B (453.01 GB)** hashed in one resumable pass;
**27/27 verified**. The human inventory of every large ignored asset — origin,
reconstructibility, backup expectation — is [docs/data-inventory.md](docs/data-inventory.md)
(571.9 GB across ten trees).

`author_log.json` answers "what is this, where did it come from, which version, when was it
taken, and is it still intact": per-source `kind`/`id`/`url`/`version`/`license`/`method`,
retrieval time, per-source file and byte counts, tree digest, declared local additions, and
the verification result and time. Sources whose origin is genuinely unknown say
`version: "unknown"` rather than guessing (`XPert-L1000`, `KPGT`, `HG`).

For the seven items whose `sources[]` has a nested asset row (CPA, chemCPA, PerturbNet,
State, TranSiGen, biolord, KPGT), the parent `source/` row's counts include that asset, so
`sources[].bytes` must not be summed to get the item total. The item total is the manifest.

---

## 3. Stage 2: four steps, two substeps

```
stage2/
├── AGENTS.md   README.md
├── 01_setup/                     nine deployed models load/train on their native assets
│   ├── AGENTS.md  README.md  model_architecture_weights_benchmarks.html
│   └── execution/<model>/{run.sh, smoke.py, env.lock, output/{*.json, *.log}}
├── 02_cm_data/
│   ├── AGENTS.md  README.md
│   ├── execution/{build_cm_datasets.py, cmdata/, tests/}
│   └── output/{corpus/v1/, datasets-1/}
├── 03_models_cmbench/            the executed 60-unit cold-drug benchmark
│   ├── AGENTS.md  README.md
│   ├── 03.1_cold_drug_18k/       primary: 18,164 genes
│   │   ├── AGENTS.md  README.md  summary.png  summary_table.csv
│   │   ├── execution/{script/, configs/, models/, tests/, Makefile}
│   │   └── output/{input/, units/, results/}
│   └── 03.2_hvg_sensitivity_5k/  sensitivity fork; same shape
├── 04_tahoe_data/                post-download processing only, no substeps
│   ├── AGENTS.md  README.md
│   ├── execution/build_pseudobulk.py
│   └── output/pseudobulk/{*.h5ad, _parts/, _state/, qa*.json, manifest.csv}
├── _envs/   audit/   _archive/   each with a README
```

- **`01_setup` deliberately has no step-level `output/`**: its artefacts are all
  intermediates and live under `execution/<model>/output/`. Recorded in ADR 0030.
- **`run_cards/` is not a directory.** Run cards are the per-model configs under
  `execution/models/<model>/`; the executed card table is
  `output/results/run_cards.csv`. Recorded in ADR 0030.
- The benchmark harness is **duplicated on purpose** in both substeps so each is
  reproducible standalone; a harness bug must be fixed in both. `repo_root.py` (duplicated)
  is the only way to find the repository root.
- Frozen benchmark evidence — `output/input/manifest.json`,
  `03.2/output/input/gene_set_5k*.json`, `output/results/{run_cards,metrics}.csv` — kept
  **byte-identical** with its historical path strings. A rename must not rewrite evidence.
- Conclusion tables (`leaderboard.csv`, `per_drug.csv`, `fidelity.csv`, `run_cards.csv`,
  `bootstrap.npz`, `comparison_18k.csv`) are force-tracked: they are the results of a
  60-unit run. The ~53 MB per-group metrics table stays local, with its path, size and
  sha256 recorded in a tracked `LOCAL_ARTIFACTS.md`.

---

## 4. The repository is now self-contained

All **285** weight/data symlinks that pointed into `/mnt/f/arc-weights` are real copies:

- 285/285 rows `ok`; **42.48 GB** added; size **and** mtime matched on every row.
- `find . -type l -lname '/mnt/*'` → **0**.
- The external archive was only read from and remains as an offline backup.
- No tracked file changed as a result — the weights stay ignored.

Before/after per-link report: `.scratch/data-provenance-restructure/ops/materialize_report.tsv`.

---

## 5. Guard rails

Three scripts under [script/](script/), all offline, self-documenting, and safe to run cold:

| Script | Question it answers | Current result |
|---|---|---|
| [check_tracked.sh](script/check_tracked.sh) | Is any executable file silently ignored by git? | PASS — 306 scanned, 335 tracked, **holes 0** |
| [check_paths.sh](script/check_paths.sh) | Does every `pert_data` / `pert_model` / `disease_data` / `drug_emb_model` path literal resolve? | PASS — 538 occurrences, 507 resolved, 31 whitelisted, **0 unresolved** |
| [verify_all.sh](script/verify_all.sh) | Does every item's `source/` still match its manifest? | PASS — **27/27**, 0 pending |

[script/path_whitelist.txt](script/path_whitelist.txt) holds 23 literal-specific exemptions,
each with a rationale: immutable ADRs, `stage2/_archive/**`, and frozen benchmark evidence.
A *live* document that merely needs updating is a path to fix, not to whitelist.

The verifiers are **fail-closed**: a missing manifest exits 2 with the command that
regenerates it, and `verify_all.sh` counts it as a failure. (Before this hardening a
deleted manifest would have reported PASS.)

---

## 6. How the migration was kept honest

Roughly 1,000 path references were rewritten. The rules that prevented silent corruption:

- **Dataset identifiers are immutable.** Only directory names and path constructions
  changed — never `temp_env_all`, `temp_env_all-smiles`, `temp_env_all_5k` as a
  `manifest.json` value, `SOURCE_DATASETS` entries or `obs` labels. Replacements were
  anchored on full path prefixes.
- **Every move is a rename, reconciled by a byte ledger.** Per item, file count and byte
  total were recorded before and after and had to match exactly, minus a stated deletion
  list. Ledgers: `ops/ticket0{3,4,5,6,7,8}-move-ledger.tsv`.
- **Nothing is deleted before its durable fact is recorded.** The manual-download note was
  archived in its own commit before being split into four items' `script/download.md`; the
  orphan dataset table became an appendix in the inventory; the model survey and validation
  verdicts were absorbed into author logs.
- **Frozen evidence keeps its historical strings.** The old-path audit is therefore over
  *tracked* files, with ADRs, archives and frozen evidence exempted — the spec's raw counts
  were all-files counts that included untracked generated artifacts.

---

## 7. Cleanup

**877.39 MB** reclaimed: throwaway recon/build workspaces (711 MB + 42 MB), 802 generated
`__pycache__` directories (114 MB), a 9.4 MB scratch tree, superseded report duplicates, a
`.bak`, caches. `.scratch/` went from 729 MB to 4.6 MB, keeping only the issue tracker and
its evidence.

---

## 8. Verification transcript

Run from the repository root; nothing here needs the network.

```
$ bash script/check_tracked.sh
scan: 306 code files on disk; 335 tracked by git
untracked: 13  (documented expected-ignored: 13, holes: 0)
RESULT: PASS — all scanned executable files are tracked (or documented expected-ignored)

$ bash script/check_paths.sh
tracked text files scanned: 753
path literals scanned: 538 occurrences, 143 distinct
resolved: 507 occurrences, 122 distinct
whitelisted: 31 occurrences, 21 distinct (entries in script/path_whitelist.txt: 23)
unresolved: 0 occurrences, 0 distinct
RESULT: PASS — every extracted path literal resolves

$ bash script/verify_all.sh
items with script/verify.sh: 27
RESULT: PASS — verified 27/27 item(s)          # ~6 min for 453 GB

$ python -m unittest discover -s stage2/02_cm_data/execution/tests
Ran 23 tests ... OK
$ python -m unittest discover -s stage2/03_models_cmbench/03.1_cold_drug_18k/execution/tests
Ran 46 tests ... OK
$ python -m unittest discover -s stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/tests
Ran 53 tests ... OK

$ find . -type l -lname '/mnt/*' | wc -l        # 0
$ find . -xtype l | wc -l                       # 0 dangling
$ git status --porcelain | wc -l                # 0
```

Every relative markdown link in the tracked docs resolves: 182 files, 169 relative links,
0 unresolved (the only exempt hits are inside the archived `stage2/_archive/PLAN.md`, which
is deliberately byte-identical).

---

## 9. Code review

A two-axis review (Standards vs. the repo's own documented contracts; Spec vs. the
originating spec) was run in parallel and everything it raised was fixed in one commit:

- **Standards: 8 findings (3 hard violations, 5 judgement calls).** Worst: five accepted
  ADRs had been silently rewritten against the index's own rule that accepted records are
  amended, never rewritten — each now carries a dated `**Update (2026-10-30, ADR 0030)**`
  line. Also fixed: ADR 0030 said each 项 has "exactly five" members when 19 of 27 have no
  `local/`; the verifier guard was fail-open; five scripts hardcoded host interpreters.
- **Spec: 4 findings (1 missing/partial, 3 wrong), no scope creep.** Worst: an `AGENTS.md`
  written by us sat directly inside `drug_emb_model/HG/source/` — repo prose inside a
  provenance `source/` tree, hashed as if it were upstream. It was moved to `HG/local/`
  and the manifest regenerated (13 → 12 files).
- Three judgement calls were **recorded rather than fixed** (they are deliberate or
  pre-existing): the duplicated root-walk helper the spec explicitly wants, `03.1`'s
  hardcoded gene count versus `03.2`'s metadata read, and the pre-existing ~100-line corpus
  build function.

Recorded in `issues/12-final-acceptance.md` → `## Code review`.

---

## 10. Divergences from the spec, as built

| Spec said | Built | Why |
|---|---|---|
| Every `scripts/` is `script/` | Our own trees yes; vendored `source/` trees keep upstream layout | `source/` must stay byte-identical to upstream |
| Paths never encode depth | Stage 2 yes; finalized Stage 1 workflows unchanged | Stage 1 was not restructured |
| `execution/{script/,configs/,run_cards/,tests/}` | `run_cards/` realised as `execution/models/<m>/` + `output/results/run_cards.csv` | the pre-restructure tree had no `run_cards/` either |
| Every Stage 2 step has `output/` | `01_setup` has none; its artefacts are intermediates under `execution/` | matches the spec's own target tree |
| Old-path counts (361 / 515 / 121) | Tracked-file audit | the spec's figures included untracked generated artifacts |

All five are written into ADR 0030's clarification section or
`ops/adjudications.md`.

---

## 11. What is deliberately not tracked, and why it still reproduces

`source/`, most of `local/`, and `stage2/**/output/` are gitignored. That is the point: the
repository tracks *how to obtain and verify* the bytes, not the bytes. A fresh clone
recovers the data by running each item's `script/download.*` and then
`script/verify_all.sh`. The tracked inventory of what is large, where it came from, whether
it can be rebuilt and how it should be backed up is
[docs/data-inventory.md](docs/data-inventory.md).

---

## 12. Next steps

1. **Merge** when ready:
   ```bash
   git checkout feat/data-provenance-restructure
   git merge --no-ff feat/data-provenance-restructure
   ```
2. **Decide on the external archive.** `/mnt/f/arc-weights` is now redundant for execution
   and kept only as an offline backup; deleting it frees ~40 GB.
3. **Decide on disk budget.** Materializing the weights moved free space on `/` from 289 GB
   to ~250 GB.
4. **If you add a new external source**, copy the shape: one 项 directory, a `script/`
   downloader, `author_log.json`, then run
   `.scratch/data-provenance-restructure/ops/hash_items.py --only <category>/<Item>` to
   produce its manifest. `script/check_tracked.sh` will tell you if you left an executable
   file untracked.
5. **Do not** hand-edit `manifest.sha256`; regenerate it. And keep `PYTHONDONTWRITEBYTECODE=1`
   in any script that imports from a `source/` tree — byte-code written into `source/` is
   contamination.

---

## 13. Where the evidence lives

| What | Where |
|---|---|
| Decision | [docs/adr/0030-external-data-provenance-and-step-layout.md](docs/adr/0030-external-data-provenance-and-step-layout.md), index in [docs/adr/README.md](docs/adr/README.md) |
| Vocabulary | [GLOSSARY.md](GLOSSARY.md) — external item/项, source tree, local addition, author log, manifest, execution area, output area, substep |
| Checks | [script/](script/) |
| Per-item provenance | `<category>/<Item>/author_log.json`, one per item; category indexes in each `README.md`/`AGENTS.md` |
| Large-asset inventory | [docs/data-inventory.md](docs/data-inventory.md) |
| Runtime evidence | `.scratch/data-provenance-restructure/ops/` — recon notes, per-ticket byte ledgers, the link-materialization report, the hashing log and tool, the docs link checker |
| Issues and resolutions | `.scratch/data-provenance-restructure/issues/01…12-*.md` (local tracker, gitignored) |
