# AGENTS.md

Operating manual for agents in this repository (AI-driven MIRI drug-repositioning, research code): how to work, where the reference docs live, and which skills to use. This is not a status board — do not add plans, progress, or result summaries here.

## Workflow

Before starting software development work in this repository, read and follow Notion’s [General Software Development Workflow](https://app.notion.com/p/413ab6347e96811b8cd9ecdaf0a972ff). It is the source of truth for the general workflow, including planning, Git, deployment, and resource monitoring. This file contains repository-specific guidance; if it conflicts with the general workflow, follow Notion.

## Execution environment (WSL2)

Before loading datasets or training, read [`docs/agents/execution-environment.md`](docs/agents/execution-environment.md) for this workstation’s shared resource limits, crash history, and data-access memory traps. Apply these repository-specific constraints alongside Notion’s resource-monitoring workflow.

## Tracker and doc conventions

- **Issue tracker**: GitHub issues at `Moke110/arc` via the `gh` CLI; specs are drafted in `.scratch/<feature-slug>/spec.md` and published. `.scratch/` is not a tracker and keeps no issue files. See `docs/agents/issue-tracker.md`.
- **Triage labels**: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.
- **Domain docs**: single-context — root `GLOSSARY.md` + `docs/adr/`. See `docs/agents/domain.md`.

## Before touching a tree

- Read the owning `AGENTS.md` (routing) and `README.md` (summary) first; every external-data category and Stage 2 step carries both.
- External data follows the four-provenance layout in ADR 0030 (`docs/adr/0030-external-data-provenance-and-step-layout.md`) — `script/`, `source/`, `local/`, plus tracked `author_log.json` and `manifest.sha256`. ADR 0030 is the law for that layout and the Stage 2 step layout. Never modify an item's `source/` tree; declared additions go in its `local/`.
- A git worktree does not share gitignored payloads: an item's `source/` and declared `local/` bytes exist only in the checkout that fetched or built them, so deleting that worktree loses them unless the main checkout already holds them (#195). Re-materialise `source/` with the item's `script/fetch.sh` then `script/verify.sh`, and `local/` with the `build_*.py` generators that `author_log.json` declares in `local_additions[]`; the tracked `local/` summaries pin the resulting digests.
- Health checks from the repo root (all offline): `script/check_paths.sh`, `script/check_tracked.sh`, `script/verify_all.sh`.

## Reference docs

- `GLOSSARY.md` — domain vocabulary; use its terms in issues, tickets, code, and docs.
- `docs/adr/` — architecture decision records; index at `docs/adr/README.md`. Read the ADRs that touch your area, and flag conflicts instead of silently overriding them.
- `docs/agents/` — issue-tracker, triage-label, domain-doc, and execution-environment rules consumed by the skills.
- `docs/data-inventory.md` — inventory of ignored large assets.
- `docs/roadmap.md` — program-level roadmap.
- `docs/wetlab-validation-stage1.md` — wet-lab validation record.
- Routing indexes: `stage1/AGENTS.md`, `stage2/01_bases/AGENTS.md`, `pert_data/AGENTS.md`, `pert_model/AGENTS.md`, `disease_data/AGENTS.md`, `drug_emb_model/AGENTS.md`; each Stage 2 step has its own `AGENTS.md` + `README.md`.
- `environment.yml` — shared conda environment `arc`; `skills-lock.json` — pinned agent-skills lockfile.

## Skills

- Project skills: `grilling`, `domain-modeling`, `code-review`, `tdd`, `diagnosing-bugs`, `prototype`, `research`, `pr`, `wizard`, `writing-for-agents`.
- Installed under `~/.agents/skills/` and symlinked from `.claude/skills/`; repo-local skills live in `.agents/skills/` and are pinned by `skills-lock.json`.
