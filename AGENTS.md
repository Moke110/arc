# AGENTS.md

Operating manual for agents in this repository (AI-driven MIRI drug-repositioning, research code): how to work, where the reference docs live, and which skills to use. This is not a status board — do not add plans, progress, or result summaries here.

## Branch rules

- **All modifications happen on `dev`.** `dev` is the working and integration branch: commit every change there. `main` only ever receives merges from `dev` — never work or commit on `main` directly.
- **Push only `main`.** `origin` receives `main` and nothing else — never push `dev` or any other branch.

## Workflow

The design-to-delivery chain is fixed; each skill consumes the previous step's output. Do not skip steps or run them out of order.

1. **`grill-with-docs`** — interrogate the plan and write decisions into the repo as they resolve: create/update the ADR under `docs/adr/`, and the affected terms in `GLOSSARY.md`.
2. **`to-spec`** — synthesize the agreed design into a spec and publish it as an issue on the tracker (`Moke110/arc`, label `ready-for-agent`).
3. **`to-tickets`** — split the spec into tracer-bullet tickets; publish one issue per ticket in dependency order with native `blocked_by` edges and `ready-for-agent`.
4. **`implement-spec`** — implement the published tickets on `dev` (the integration branch), working the frontier of the ticket graph, and finish with the `code-review` pass.

Execution rules:

- **Subagent context budget: under ~20% per agent.** Fan work out to subagents — one ticket (or one narrow prefactor scope) per implementer subagent, plus separate exploration and merger subagents. Pass context pointers (spec path, issue number, research-note path, commit SHA), not pasted content; when a window fills, hand off to a fresh subagent.
- **Close the tracker at the end: no open issue survives the run.** Resolve every ticket: completed work closes with its evidence; unfinished work is recorded as a comment on the owning ticket (what remains and why) and the ticket is then closed. Report the `dev` state and the closed-ticket list.

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
- `docs/agents/` — issue-tracker, triage-label, and domain-doc rules consumed by the skills.
- `docs/data-inventory.md` — inventory of ignored large assets.
- `docs/roadmap.md` — program-level roadmap.
- `docs/wetlab-validation-stage1.md` — wet-lab validation record.
- Routing indexes: `stage1/AGENTS.md`, `stage2/AGENTS.md`, `pert_data/AGENTS.md`, `pert_model/AGENTS.md`, `disease_data/AGENTS.md`, `drug_emb_model/AGENTS.md`; each Stage 2 step has its own `AGENTS.md` + `README.md`.
- `environment.yml` — shared conda environment `arc`; `skills-lock.json` — pinned agent-skills lockfile.

## Skills

- Workflow chain: `grill-with-docs` (grilling + domain-modeling), `to-spec`, `to-tickets`, `implement-spec`.
- Supporting: `grilling`, `domain-modeling`, `code-review`, `tdd`, `diagnosing-bugs`, `prototype`, `research`, `pr`, `wizard`, `writing-for-agents`.
- Installed under `~/.agents/skills/` and symlinked from `.claude/skills/`; repo-local skills live in `.agents/skills/` and are pinned by `skills-lock.json`.
