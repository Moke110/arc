"""
Entrant audit (ticket 02, re-run) — a survey of the Stage 2 candidate models.

Re-audit criteria (these supersede the pretrained-weights-only rule of the original
ticket-02 audit):

- **Capability gates (hard):** train/predict on the full 20,570-gene universe, predict a
  cold drug from chemical structure, and take a dose for the cold-dose column.
- **Compute tiers:** Tier L (>50 M params, or a pretraining corpus not reproducible here)
  must ship a public checkpoint. Tier S may train locally from scratch on the fold's train
  wells, with a total benchmark-training budget of 16 h on the RTX 4070 Ti SUPER.
- **License gate:** permissive or explicit research permission; otherwise excluded.

Each record carries the evidence used at audit time (repo commit, license, checkpoint,
native gene space, cold-drug/dose capability, estimated hours) and an explicit verdict:
`admitted`, `reference` (runnable but not a contender), or `excluded`.
"""

import csv
from pathlib import Path

# tier: "S" (local training allowed) | "L" (public checkpoint required)
# status: "admitted" | "reference" | "excluded"
AUDIT = [
    # ------------------------------------------------------------------ admitted
    {
        "model": "chemCPA",
        "tier": "S",
        "params": "est. 1-5M",
        "repo": "https://github.com/theislab/chemCPA",
        "license": "MIT",
        "pinned_commit": "43e830eb0958c54e4aa64442c17ec0fed19b3f15",
        "checkpoint_source": "https://f003.backblazeb2.com/file/chemCPA-models/chemCPA_models.zip",
        "checkpoint_form": "raw ComPert state dicts in the public zip: 977 or 2,000 genes, "
                           "188 drugs, 3 covariates; no 5,000- or 20,570-gene checkpoint. "
                           "The drug encoder is gene-agnostic, but the decoder is not.",
        "gene_space": "977 / 2,000 (public checkpoints); 20,570 only after retraining",
        "cold_drug": "yes",
        "dose": "yes",
        "est_hours": 10,
        "recipe": "public checkpoint proves load/forward on L1000; for a 20,570-gene CM run, "
                  "reuse the chemical encoder and retrain the VAE decoder on the fold's CM "
                  "train wells with drug = chemical embedding and condition = per-plate "
                  "baseline; predict held-out drugs",
        "status": "admitted",
        "reason": "chemical encoder generalises to unseen drugs and the model is dose-aware; "
                  "the public checkpoints are L1000-only, so whole-transcriptome use requires "
                  "Tier S retraining inside the 16 h budget",
    },
    {
        "model": "CPA",
        "tier": "S",
        "params": "est. ~1M",
        "repo": "https://github.com/theislab/CPA",
        "license": "BSD-3-Clause",
        "pinned_commit": "fbd7c0250edc23eff003a10c99655579c53afd63",
        "checkpoint_source": "Google Drive folder linked from the CPA README "
                             "('Datasets and Pre-trained models')",
        "checkpoint_form": "models linked from the README; not verified for reuse at the "
                           "20,570-gene universe, so this entrant is trained locally",
        "gene_space": "20,570 (VAE gene dimension)",
        "cold_drug": "yes",
        "dose": "yes",
        "est_hours": 9,
        "recipe": "train CPA on each fold's CM train wells; reuse chemCPA's chemical encoder "
                  "for the held-out drug (the CPA paper's unseen-drug result uses that "
                  "extension); dose is a native covariate",
        "status": "admitted",
        "reason": "CPA's OOD result covers completely unseen drugs via a chemical encoder; "
                  "per-dataset training is minutes-to-hours on the target GPU",
    },
    {
        "model": "biolord",
        "tier": "S",
        "params": "est. 1-5M",
        "repo": "https://github.com/nitzanlab/biolord",
        "license": "BSD-3-Clause",
        "pinned_commit": "b7688790e49728d7f8d0906b980a629de484b19b",
        "checkpoint_source": "",
        "checkpoint_form": "none found in the repo/PyPI docs at audit time; per-dataset training only",
        "gene_space": "20,570 (VAE gene dimension)",
        "cold_drug": "yes",
        "dose": "yes (ordered attribute)",
        "est_hours": 12,
        "recipe": "train biolord on each fold's CM train wells; drug = RDKit chemically "
                  "informed embedding, dose = ordered attribute, baseline well = control "
                  "cell; counterfactual prediction for held-out drugs",
        "status": "admitted",
        "reason": "the paper's OOD split holds out nine unseen drugs and predicts them from "
                  "chemical features; single VAE, tractable on one 16 GB card",
    },
    {
        "model": "PRnet",
        "tier": "S",
        "params": "~1-3M (checkpoints 2.2 MB / 8.4 MB)",
        "repo": "https://github.com/Perturbation-Response-Prediction/PRnet",
        "license": "Apache-2.0",
        "pinned_commit": "f19174bde3ed2633f54c7831799cc38c4ffc7a0d",
        "checkpoint_source": "in-repo checkpoint/lincs_best_epoch_all.pt and "
                             "sciplex_best_epoch_all.pt; also Zenodo 10.5281/zenodo.13751384",
        "checkpoint_form": "published checkpoints are 978-landmark (bulk) / 5,000-HVG "
                           "(single-cell); not reusable for the 20,570-gene universe",
        "gene_space": "20,570 (x_dimension is configurable; the published config is 978/5,000)",
        "cold_drug": "yes",
        "dose": "yes",
        "est_hours": 8,
        "recipe": "set x_dimension=20,570 in the training config and train on each fold's CM "
                  "train wells; drug = rFCFP chemical embedding, dose = native dose encoder; "
                  "predict held-out drugs",
        "status": "admitted",
        "reason": "explicitly designed for novel compounds and dosages; Apache-2.0; tiny "
                  "model, so a 20,570-gene re-train fits the budget",
    },
    # ------------------------------------------------------------------ reference
    {
        "model": "State",
        "tier": "L",
        "params": "600M (SE-600M) + transition model",
        "repo": "https://github.com/ArcInstitute/state",
        "license": "Arc Research Institute State Model Non-Commercial License",
        "pinned_commit": "9bbfe78a434a55205e4de834e1ea99f85f7a3add",
        "checkpoint_source": "Arc Institute State pretrained weights (SE-600M + ST-HVG-Tahoe)",
        "checkpoint_form": "cell-embedding *.ckpt + transition run dir",
        "gene_space": "2,000 HVG (1,398 in the common universe)",
        "cold_drug": "no",
        "dose": "yes (inside the Tahoe one-hot)",
        "est_hours": 16,
        "recipe": "already implemented in 03_state; keep the cold-dose row only",
        "status": "reference",
        "reason": "perturbation identity is a fixed 1,138-way Tahoe drug-dose one-hot: 62/79 "
                  "CM drugs are unseen and fall back to control, and there is no CM cell type. "
                  "Kept as a reference row, not a contender.",
    },
    # ------------------------------------------------------------------ excluded
    {
        "model": "scGen",
        "tier": "S",
        "params": "est. 1-2M",
        "repo": "https://github.com/theislab/scgen",
        "license": "GPL-3.0",
        "pinned_commit": "d79e1f04233c30f9a4eb5b8d57718127909807d7",
        "checkpoint_source": "",
        "checkpoint_form": "none found in the repo/PyPI docs at audit time",
        "gene_space": "20,570 (VAE gene dimension)",
        "cold_drug": "no",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "scGen transfers a latent perturbation shift to unseen cell types; it has "
                  "no drug encoder and no dose input, so it cannot predict a held-out drug",
    },
    {
        "model": "CellOT",
        "tier": "S",
        "params": "est. ~1M (ICNN; not reported)",
        "repo": "https://github.com/bunnech/cellot",
        "license": "BSD-3-Clause",
        "pinned_commit": "522d2b953da8ad244fcf36f64521487fdf763788",
        "checkpoint_source": "",
        "checkpoint_form": "none; one OT map is fitted per perturbation",
        "gene_space": "20,570 (input dimension)",
        "cold_drug": "no",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "CellOT learns a separate transport map per perturbation and generalises "
                  "to unseen patients/cells, not unseen drugs; it has no chemical or dose "
                  "input",
    },
    {
        "model": "scPRAM",
        "tier": "S",
        "params": "est. 1-10M (VAE + attention; not reported)",
        "repo": "https://github.com/jiang-q19/scPRAM",
        "license": "MIT",
        "pinned_commit": "3f490343d9703ba4b7223ec088b0b5a24338255e",
        "checkpoint_source": "",
        "checkpoint_form": "checkpointing during training only; no published pretrained weights",
        "gene_space": "20,570 (VAE gene dimension)",
        "cold_drug": "no",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "scPRAM's stated generalisation is to unseen cell types; the perturbation "
                  "is a categorical embedding with no chemical structure or dose input",
    },
    {
        "model": "CFM-GP",
        "tier": "S",
        "params": "est. <1M (MLP vector field)",
        "repo": "https://github.com/abrarrahmanabir/CFM-GP",
        "license": "GPL-3.0",
        "pinned_commit": "3046375ba645cb4cfbdebfbfb375ca879f648abd",
        "checkpoint_source": "README lists pretrained models",
        "checkpoint_form": "per-cell-type flow-matching checkpoints",
        "gene_space": "20,570 (input dimension)",
        "cold_drug": "no",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "the vector field is conditioned on cell type only (gene perturbation "
                  "across cell types); there is no drug identity or dose in the model",
    },
    {
        "model": "TxPert",
        "tier": "L",
        "params": "not reported (graph transformer over knowledge graphs)",
        "repo": "https://github.com/valence-labs/TxPert",
        "license": "Recursion Non-Commercial End User License (CC BY-NC-SA-like)",
        "pinned_commit": "08d82eea86746b044cf7531f4ec8c5f60e1cb73f",
        "checkpoint_source": "in-repo weights for the top public-data variants",
        "checkpoint_form": "latent-transfer checkpoint over gene knowledge graphs",
        "gene_space": "full transcriptome (gene nodes)",
        "cold_drug": "no",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "TxPert predicts genetic perturbations (CRISPR KO/OE); a drug would first "
                  "have to be mapped to target genes, which changes the task and is not a "
                  "chemical-structure cold-drug prediction",
    },
    {
        "model": "PertDiff (roadmap: MCDP)",
        "tier": "S",
        "params": "not reported (diffusion MLP + graph encoder + LLM embeddings)",
        "repo": "https://github.com/Panda-myj/PertDiff",
        "license": "none declared in the repo (Zenodo license not verified)",
        "pinned_commit": "073a92dd154bf096cdb8878f60ffa3add59ae75d",
        "checkpoint_source": "Zenodo 10.5281/zenodo.18427848 (model.pt)",
        "checkpoint_form": "trained diffusion checkpoint",
        "gene_space": "978 landmark genes (then a linear expansion to 10,174)",
        "cold_drug": "yes",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "the model is 978-landmark by construction (attention shapes are hard-coded "
                  "at 978) and the repo declares no license; both violate the Stage 2 gates",
    },
    {
        "model": "LDM-Drug",
        "tier": "S",
        "params": "not reported (VAE + latent diffusion)",
        "repo": "https://github.com/bmil-jnu/Perturbed-GE-LDM",
        "license": "Pre-publication Research License (view-and-study only; no use without "
                   "written permission)",
        "pinned_commit": "7814c21437b9ecfe82aa21c68ad4e3da5a622a04",
        "checkpoint_source": "Google Drive (best_ldm.pt), linked from the README",
        "checkpoint_form": "VAE + diffusion checkpoint",
        "gene_space": "978 landmark genes only",
        "cold_drug": "yes",
        "dose": "yes",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "trained and evaluated on 978 L1000 landmarks only, and the license forbids "
                  "use without written permission; excluded until a 20,570-gene refactor and "
                  "a usable license exist",
    },
    {
        "model": "scGPT",
        "tier": "L",
        "params": "~33-53M",
        "repo": "https://github.com/bowang-lab/scGPT",
        "license": "MIT",
        "pinned_commit": "cebd6fae655b9c585a4807daa3ac31bb764f06b4",
        "checkpoint_source": "public scGPT checkpoints (figshare/Hugging Face)",
        "checkpoint_form": "whole-human / organ pretrained transformer",
        "gene_space": "1,200 input HVGs in the perturbation tasks (paper)",
        "cold_drug": "no",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "the perturbation head is gene-token based with a ~1,200-HVG input cap and "
                  "no chemical or dose conditioning; adding a drug encoder would be a new "
                  "model, not a fine-tune",
    },
    {
        "model": "scFoundation",
        "tier": "L",
        "params": "100M",
        "repo": "https://github.com/biomap-research/scFoundation",
        "license": "Apache-2.0",
        "pinned_commit": "397631c495eddf9ad6644fc00c6ea8139e651245",
        "checkpoint_source": "public scFoundation checkpoint (~19k genes)",
        "checkpoint_form": "cell-embedding foundation model",
        "gene_space": "~19,000 genes (close to the universe)",
        "cold_drug": "no",
        "dose": "no",
        "est_hours": 0,
        "recipe": "",
        "status": "excluded",
        "reason": "a cell-embedding foundation model with no drug-conditioned perturbation "
                  "head; using it would mean training a new drug encoder and decoder head, "
                  "not fine-tuning the released model",
    },
]


def admitted():
    return [r for r in AUDIT if r["status"] == "admitted"]


def reference():
    return [r for r in AUDIT if r["status"] == "reference"]


def excluded():
    return [r for r in AUDIT if r["status"] == "excluded"]


def check_audit_gate():
    """Fail loudly on a verdict that violates the re-audit criteria."""
    for r in AUDIT:
        assert r["status"] in ("admitted", "reference", "excluded"), r["model"]
        if r["status"] == "admitted":
            if r["cold_drug"] != "yes":
                raise AssertionError(f"{r['model']}: admitted without cold-drug capability")
            if "20,570" not in r["gene_space"]:
                raise AssertionError(
                    f"{r['model']}: admitted without a 20,570-gene output space")
            if r["est_hours"] > 16:
                raise AssertionError(f"{r['model']}: estimate exceeds the 16 h budget")
            if not r["license"] or r["license"].lower() in ("none", "unclear"):
                raise AssertionError(f"{r['model']}: admitted without a usable license")
            if r["tier"] == "L" and not r["checkpoint_source"]:
                raise AssertionError(
                    f"{r['model']}: Tier L admitted without a public checkpoint")
            if r["tier"] == "S" and not r["recipe"]:
                raise AssertionError(
                    f"{r['model']}: Tier S admitted without a local-training recipe")
        else:
            if not r["reason"]:
                raise AssertionError(f"{r['model']}: {r['status']} without a reason")
    return True


def write_audit(out_csv, out_md):
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(AUDIT[0].keys()))
        w.writeheader()
        w.writerows(AUDIT)

    lines = [
        "# Stage 2 entrant audit — candidate models (ticket 02 re-run)",
        "",
        "Capability gates come first: a whole-transcriptome output on the 20,570-gene "
        "universe, a cold-drug prediction from chemical structure, and a dose input for the "
        "cold-dose column. **Tier S** models may be trained locally on the fold's train wells "
        "(≤16 h for the full benchmark on the RTX 4070 Ti SUPER); **Tier L** models require a "
        "public pretrained checkpoint. Evidence was read from the cited papers and repos at "
        "audit time. This survey supersedes the pre-re-audit status lines in "
        "`03_*/README.md`; the verdicts below are the model-selection record.",
        "",
        "| model | tier | status | cold-drug | dose | gene space | est h | license | checkpoint |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in AUDIT:
        lines.append(
            f"| {r['model']} | {r['tier']} | {r['status']} | {r['cold_drug']} | {r['dose']} | "
            f"{r['gene_space']} | {r['est_hours'] or '—'} | {r['license']} | "
            f"{r['checkpoint_source'] or '—'} |")

    lines += ["", "## Admitted (usable entrants)", ""]
    for r in admitted():
        lines.append(
            f"- **{r['model']}** ({r['tier']} tier, ~{r['est_hours']} h) — {r['reason']}. "
            f"Recipe: {r['recipe']}.")
    lines += ["", "## Reference only (runnable, not a contender)", ""]
    for r in reference():
        lines.append(f"- **{r['model']}** — {r['reason']}")
    lines += ["", "## Excluded", ""]
    for r in excluded():
        lines.append(f"- **{r['model']}** — {r['reason']}")

    lines += ["", "## Shared adapter risks", "",
              "- **Input scale.** The benchmark freezes log2(CPM+1), while chemCPA / CPA / "
              "biolord natively model counts. Each adapter must either invert to counts or "
              "switch to a Gaussian likelihood; the fold-0 pilot must show a sane "
              "reconstruction before the full run.",
              "- **Bulk-as-cell.** DToxS and DToxS-iPSC are bulk 3'-DGE and gse244740 is "
              "pseudo-bulk. The single-cell VAEs treat each well as one cell; this is an "
              "approximation to document, not a hidden one.",
              "- **Drug identity.** chemCPA / CPA need a chemical encoder, biolord needs "
              "RDKit features, PRnet uses rFCFP. All consume SMILES, but the encoder must be "
              "pinned and its gene-independence verified.",
              "- **RAM, not VRAM.** The host has 15 GB RAM and 16 GB VRAM. The precomputed "
              "bench matrices are ~0.4 GB, so training fits; raw h5ad loading must be "
              "streamed/chunked."]

    lines += ["", "## What would change a verdict", "",
              "- A State adapter that re-initialises the perturbation embedding could improve "
              "cold-dose, but cold-drug remains impossible without chemical initialisation.",
              "- LDM-Drug / PertDiff could re-enter after a 20,570-gene refactor and a usable "
              "license; that is a new model, not the published checkpoint.",
              "- scGPT / scFoundation could re-enter only with a drug-conditioned head trained "
              "on CM data — a research project, not a fine-tune.",
              ""]
    out_md = Path(out_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    while lines and lines[-1] == "":
        lines.pop()
    out_md.write_text("\n".join(lines) + "\n")
