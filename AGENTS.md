# CarPS agent instructions

## Start work

1. At the start of each session, use the Notion connector to read the full [Agent通用工作流 General agent workflow](https://app.notion.com/p/413ab6347e96811b8cd9ecdaf0a972ff) before giving a substantive reply or executing the task. Follow it for planning, Git, deployment and resource monitoring. Use these repository rules alongside it and resolve general workflow conflicts in favor of the Notion page. If the connector is unavailable, denies access or returns incomplete text, briefly say so in your first reply and continue with the available instructions.
2. Run repo-relative commands from the root of the checkout you are working in. Keep disposable files and working directories under that checkout's ignored `tmp/`. Point tools that use `TMPDIR` there. Keep issue-spec drafts at `.scratch/<feature-slug>/spec.md` under the separate tracker convention below.
3. Before loading data, training or running heavy verification, read the [WSL2 execution guide](docs/agents/execution-environment.md). Check live WSL and Windows host resources, cap heavy jobs and account for concurrent worktrees. Use the `arc` conda environment where the owning guide calls for it. [Environment definition](environment.yml)
4. Before changing a tree, read its owning `AGENTS.md` and `README.md`. Use the [glossary](GLOSSARY.md) for names and terms in code, docs and issues. Read the relevant [ADRs](docs/adr/README.md) and flag conflicts with accepted decisions. [Domain document conventions](docs/agents/domain.md)

## Route by owner

| Work | Owning guide |
| --- | --- |
| CarPS-Transcriptome (CPST) data, transcriptome model or experiments | [CPST guide](cpst/AGENTS.md) and [report](cpst/README.md) |
| CarPS-Phenotype (CPSP) corpus, phenotype model or benchmark | [CPSP guide](cpsp/AGENTS.md) and [report](cpsp/README.md) |
| Shared CPST training, evaluation and plotting | [Utility guide](cpst/util/AGENTS.md) and [API inventory](cpst/util/README.md) |
| External perturbation datasets | [pert_data guide](pert_data/AGENTS.md) |
| External perturbation models | [pert_model guide](pert_model/AGENTS.md) |
| External drug embedding libraries | [drug_emb_model guide](drug_emb_model/AGENTS.md) |
| External disease transcriptomes | [disease_data guide](disease_data/AGENTS.md) |
| External phenotype assays | [pheno_data guide](pheno_data/AGENTS.md) |

The Stage 1 predecessor in `stage1/` is private and local. Consult its local guide for work on the predecessor. Keep public workflows and required checks runnable with that tree absent. [Public release boundary](docs/adr/0106-cpst-and-cpsp-public-release-boundaries.md)

## External data and verification

- Preserve upstream bytes in `source/` and put declared additions in `local/`. Use the item's `script/`, tracked `author_log.json` and `manifest.sha256` for acquisition and verification. [External provenance contract](docs/adr/0030-external-data-provenance-and-step-layout.md)
- Each worktree has its own ignored payloads. Before removing a worktree, check where needed inputs and generated products will remain available. Materialize missing inputs in the checkout you will use. Follow the owning acquisition instructions, then run the item's `script/verify.sh`. Rebuild local additions with the generators declared in `author_log.json` under `local_additions[]` and check the tracked digests where supplied. [Data inventory](docs/data-inventory.md)
- Run the applicable offline checks with Bash from the checkout root. Use the local payloads and the WSL2 limits for the full integrity check. [External path audit](script/check_paths.sh) · [Tracked-code audit](script/check_tracked.sh) · [External integrity checks](script/verify_all.sh) · [CPST path audit](script/check_cpst_paths.sh) · [CPSP path audit](script/check_cpsp_paths.sh)

## Tracker and skills

Use `gh` for issues in `Moke110/carps`. Publish spec drafts there and keep local drafts in the spec workspace. Apply the project labels when triaging. [Tracker conventions](docs/agents/issue-tracker.md) · [Triage labels](docs/agents/triage-labels.md)

Read the applicable `SKILL.md` before using a project skill. Shared skills live under `~/.agents/skills/`, and repository skills under `.agents/skills/`. [Skill lockfile](skills-lock.json)
