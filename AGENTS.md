# AGENTS.md

Operating manual for agents in this repository (AI-driven MIRI drug-repositioning, research code): how to work, where the reference docs live, and which skills to use. This is not a status board — do not add plans, progress, or result summaries here.

## Branch rules

- **All modifications happen on `dev`.** `dev` is the working and integration branch: commit every change there. `main` only ever receives merges from `dev` — never work or commit on `main` directly.
- **Push only `main`.** `origin` receives `main` and nothing else — never push `dev` or any other branch.

## Execution environment (WSL2)

This repo is developed inside one WSL2 distro on a single workstation, so every agent,
worktree and subagent shares the same memory. Measured envelope and the crash record:
[`docs/agents/execution-environment.md`](docs/agents/execution-environment.md). Read it
before any job that loads datasets or trains.

- **Measured now** (48 GB `.wslconfig` active since the 14:35 restart): 8 vCPU, **47.0 GiB
  RAM**, 64 GiB swap on `F:`, one 16 GB GPU, 207 GB free disk. Budget against `free -h`, not
  against `.wslconfig` — a config change needs a `wsl --shutdown`.
- **Two heavy jobs repo-wide, no more.** Count worktrees and subagents together; keep their
  combined peak under ~28 GiB. Six concurrent implementers is how the 2026-11-02 run crashed
  this VM (the pre-restart 31 GiB envelope). Long training runs go one at a time.
- **Cap each heavy process**: `systemd-run --user --scope -p MemoryMax=10G
  -p MemorySwapMax=8G -- <cmd>` (18G for a single formal pretrain, 12G for the 06.3/06.4
  data jobs), and check `free -g` first — abort if available is under 12 GiB.
- **Bound the memory of data access**: no dense observation × 18,164 matrix, no
  `anndata.read_h5ad(backed="r")` on a large file (it materialises `layers`: 8.1-12.9 GB —
  see the trap in the environment doc), `h5py` header reads + sliced/CSR access, `float32`,
  chunk per fold/dose, `.npy` memmap for large caches, never `n_jobs=-1`; with two jobs
  running set `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=2`.
- **Record peaks**: run >1 GiB jobs under `/usr/bin/time -v` and put the peak RSS in the
  ticket record; a rising peak in a test suite or build gate is a defect.

## Workflow

Before starting software development work in this repository, read and follow Notion’s [软件通用工作流](https://app.notion.com/p/413ab6347e96811b8cd9ecdaf0a972ff).

## Tracker and doc conventions

- **Issue tracker**: GitHub issues at `Moke110/arc` via the `gh` CLI; specs are drafted in `.scratch/<feature-slug>/spec.md` and published. `.scratch/` is not a tracker and keeps no issue files. See `docs/agents/issue-tracker.md`.
- **Triage labels**: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.
- **Domain docs**: single-context — root `GLOSSARY.md` + `docs/adr/`. See `docs/agents/domain.md`.

## Before touching a tree

- Read the owning `AGENTS.md` (routing) and `README.md` (summary) first; every external-data category and Stage 2 step carries both.
- External data follows the four-provenance layout in ADR 0030 (`docs/adr/0030-external-data-provenance-and-step-layout.md`) — `script/`, `source/`, `local/`, plus tracked `author_log.json` and `manifest.sha256`. ADR 0030 is the law for that layout and the Stage 2 step layout. Never modify an item's `source/` tree; declared additions go in its `local/`.
- Health checks from the repo root (all offline): `script/check_paths.sh`, `script/check_tracked.sh`, `script/verify_all.sh`.

## Reference docs

- `GLOSSARY.md` — domain vocabulary; use its terms in issues, tickets, code, and docs.
- `docs/adr/` — architecture decision records; index at `docs/adr/README.md`. Read the ADRs that touch your area, and flag conflicts instead of silently overriding them.
- `docs/agents/` — issue-tracker, triage-label, domain-doc, and execution-environment rules consumed by the skills.
- `docs/data-inventory.md` — inventory of ignored large assets.
- `docs/roadmap.md` — program-level roadmap.
- `docs/wetlab-validation-stage1.md` — wet-lab validation record.
- Routing indexes: `stage1/AGENTS.md`, `stage2/AGENTS.md`, `pert_data/AGENTS.md`, `pert_model/AGENTS.md`, `disease_data/AGENTS.md`, `drug_emb_model/AGENTS.md`; each Stage 2 step has its own `AGENTS.md` + `README.md`.
- `environment.yml` — shared conda environment `arc`; `skills-lock.json` — pinned agent-skills lockfile.

## Skills

- Project skills: `grilling`, `domain-modeling`, `code-review`, `tdd`, `diagnosing-bugs`, `prototype`, `research`, `pr`, `wizard`, `writing-for-agents`.
- Installed under `~/.agents/skills/` and symlinked from `.claude/skills/`; repo-local skills live in `.agents/skills/` and are pinned by `skills-lock.json`.
