"""
Checkpoint-availability audit (ticket 02) — pre-flight entrant gate.

Stage 2 v1 is compute-constrained: every entrant must start from a **public
pretrained checkpoint** and be adapted to CM data only by **finetuning**.
Nothing is pretrained or trained from scratch. This module records, per
candidate, its public repository (pinned commit + license), the public
checkpoint and the form it ships in, and a CM-finetune recipe — or an explicit
reason the candidate is excluded.

The `gate` fails if any admitted entrant has no checkpoint or no finetune
recipe, and if any candidate is neither admitted nor explicitly excluded.

Repository/commit/license facts were recorded from the public repos at audit
time; exact checkpoint hosting and recipes are implementation details (the spec
fixes the policy, not the URLs).
"""

import csv
from pathlib import Path

# status: "admitted" | "excluded"
AUDIT = [
    {
        "model": "chemCPA",
        "repo": "https://github.com/theislab/chemCPA",
        "license": "MIT",
        "pinned_commit": "43e830eb0958c54e4aa64442c17ec0fed19b3f15",
        "checkpoint_source": "https://f003.backblazeb2.com/file/chemCPA-models/chemCPA_models.zip",
        "checkpoint_form": "zip of trained chemCPA models (README 'weight checkpoints')",
        "finetune_recipe": "load a public checkpoint, then train on the fold's CM "
                           "train wells (drug = chemical embedding, condition = baseline "
                           "transcriptome adapter); predict held-out drugs",
        "status": "admitted",
        "reason": "public pretrained weights + documented training loop",
    },
    {
        "model": "State",
        "repo": "https://github.com/ArcInstitute/state",
        "license": "Arc Research Institute State Model Non-Commercial License",
        "pinned_commit": "9bbfe78a434a55205e4de834e1ea99f85f7a3add",
        "checkpoint_source": "Arc Institute State pretrained weights (SE-600M cell "
                             "model + State transition model), via the `state` CLI",
        "checkpoint_form": "cell-embedding checkpoint (*.ckpt) + transition run dir "
                           "(config.yaml + checkpoints/)",
        "finetune_recipe": "`state emb` (or load the SE-600M checkpoint) to embed CM "
                           "baselines, then `state tx train` on the fold's CM train "
                           "wells; `state tx predict` for held-out drugs",
        "status": "admitted",
        "reason": "public pretrained cell + transition checkpoints with a finetuning path",
    },
    {
        "model": "CPA",
        "repo": "https://github.com/theislab/CPA",
        "license": "BSD-3-Clause",
        "pinned_commit": "fbd7c0250edc23eff003a10c99655579c53afd63",
        "checkpoint_source": "",
        "checkpoint_form": "",
        "finetune_recipe": "",
        "status": "excluded",
        "reason": "no public pretrained checkpoint published in the repo/docs; CPA is "
                  "trained per dataset, i.e. from scratch — violates the compute "
                  "constraint",
    },
    {
        "model": "biolord",
        "repo": "https://github.com/nitzanlab/biolord",
        "license": "BSD-3-Clause",
        "pinned_commit": "b7688790e49728d7f8d0906b980a629de484b19b",
        "checkpoint_source": "",
        "checkpoint_form": "",
        "finetune_recipe": "",
        "status": "excluded",
        "reason": "no public pretrained checkpoint; biolord is trained per dataset, "
                  "i.e. from scratch — violates the compute constraint",
    },
    {
        "model": "scGen",
        "repo": "https://github.com/theislab/scgen",
        "license": "GPL-3.0",
        "pinned_commit": "d79e1f04233c30f9a4eb5b8d57718127909807d7",
        "checkpoint_source": "",
        "checkpoint_form": "",
        "finetune_recipe": "",
        "status": "excluded",
        "reason": "no public pretrained checkpoint in the repo; scGen is trained per "
                  "dataset, i.e. from scratch — violates the compute constraint",
    },
]


def admitted():
    return [r for r in AUDIT if r["status"] == "admitted"]


def excluded():
    return [r for r in AUDIT if r["status"] == "excluded"]


def check_audit_gate():
    """Fail loudly on an entrant with neither a checkpoint nor a reason to drop."""
    for r in AUDIT:
        assert r["status"] in ("admitted", "excluded"), r["model"]
        if r["status"] == "admitted":
            if not r["checkpoint_source"]:
                raise AssertionError(
                    f"{r['model']}: admitted without a public checkpoint")
            if not r["finetune_recipe"]:
                raise AssertionError(
                    f"{r['model']}: admitted without a CM-finetune recipe")
        else:
            if not r["reason"]:
                raise AssertionError(f"{r['model']}: excluded without a reason")
    return True


def write_audit(out_csv, out_md):
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(AUDIT[0].keys()))
        w.writeheader()
        w.writerows(AUDIT)
    lines = ["# Stage 2 v1 entrant audit (ticket 02)", "",
             "Pre-flight checkpoint-availability gate. Admitted entrants have a "
             "public pretrained checkpoint and a CM-finetune recipe; excluded "
             "candidates are recorded with reasons. Stage 2 v1 is "
             "pretrained-weights-only (no from-scratch training).", "",
             "| model | status | repo | license | pinned commit | checkpoint |",
             "|---|---|---|---|---|---|"]
    for r in AUDIT:
        lines.append(
            f"| {r['model']} | {r['status']} | [repo]({r['repo']}) | "
            f"{r['license']} | `{r['pinned_commit'][:10]}` | "
            f"{r['checkpoint_source'] or '—'} |")
    lines += ["", "## Exclusions", ""]
    for r in excluded():
        lines.append(f"- **{r['model']}** — {r['reason']}")
    out_md = Path(out_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n")
