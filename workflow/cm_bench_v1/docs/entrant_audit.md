# Stage 2 v1 entrant audit (ticket 02)

Pre-flight checkpoint-availability gate. Admitted entrants have a public pretrained checkpoint and a CM-finetune recipe; excluded candidates are recorded with reasons. Stage 2 v1 is pretrained-weights-only (no from-scratch training).

| model | status | repo | license | pinned commit | checkpoint |
|---|---|---|---|---|---|
| chemCPA | admitted | [repo](https://github.com/theislab/chemCPA) | MIT | `43e830eb09` | https://f003.backblazeb2.com/file/chemCPA-models/chemCPA_models.zip |
| State | admitted | [repo](https://github.com/ArcInstitute/state) | Arc Research Institute State Model Non-Commercial License | `9bbfe78a43` | Arc Institute State pretrained weights (SE-600M cell model + State transition model), via the `state` CLI |
| CPA | excluded | [repo](https://github.com/theislab/CPA) | BSD-3-Clause | `fbd7c0250e` | — |
| biolord | excluded | [repo](https://github.com/nitzanlab/biolord) | BSD-3-Clause | `b7688790e4` | — |
| scGen | excluded | [repo](https://github.com/theislab/scgen) | GPL-3.0 | `d79e1f0423` | — |

## Exclusions

- **CPA** — no public pretrained checkpoint published in the repo/docs; CPA is trained per dataset, i.e. from scratch — violates the compute constraint
- **biolord** — no public pretrained checkpoint; biolord is trained per dataset, i.e. from scratch — violates the compute constraint
- **scGen** — no public pretrained checkpoint in the repo; scGen is trained per dataset, i.e. from scratch — violates the compute constraint
