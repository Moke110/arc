"""Raw-first drug representation comparison connected to the shared Stage 2 releases."""

import json
from pathlib import Path
import sys
from tempfile import NamedTemporaryFile, TemporaryDirectory

import numpy as np
import pandas as pd

STAGE2 = Path(__file__).resolve().parents[4]
DOSE_API = str(STAGE2 / "07_tsg_dose/output/api")
if DOSE_API not in sys.path:
    sys.path.insert(0, DOSE_API)

from tsg_raw_first.evaluation import score
from tsg_raw_first.inputs import Inputs
from stage2_util.representations import ARMS, RepresentationTable, WIDTHS
from stage2_util.seeds import read_seeds

SCHEMA = "tsg-drug-emb-raw-first-v1"
SCENARIO = "cold-drug"


def _bind(inputs, representation, representation_sha256):
    data = Inputs(inputs)
    table = RepresentationTable(representation, source=data.source,
                                manifest_sha256=representation_sha256)
    if table.manifest_sha256 != data.config["representation_sha256"]:
        raise ValueError("representation generation differs from the paired input")
    table.verified_parent_keys(data.parents)
    seed = data.config["training_seed"]
    if seed not in read_seeds():
        raise ValueError("representation comparison training seed differs from the paired input")
    provenance = {"schema": SCHEMA, "source": data.source, "scenario": SCENARIO,
                  "pair_sha256": data.pair_sha256, "view_generation": data.config["view_generation"],
                  "view_payload_sha256": data.config["view_payload_sha256"],
                  "representation_sha256": table.manifest_sha256, "representations": WIDTHS,
                  "seed": seed}
    return data, table, provenance


def prepare(inputs, receipt, *, representation=None, representation_sha256=None):
    """Verify one paired source and bind its parent identities to all shared representations."""
    _, _, provenance = _bind(inputs, representation, representation_sha256)
    receipt = Path(receipt)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", dir=receipt.parent, encoding="utf-8") as staged:
        json.dump(provenance, staged, indent=2, sort_keys=True)
        staged.write("\n")
        staged.flush()
        receipt.hardlink_to(staged.name)
    return provenance


class Comparison:
    """Open the exact raw-first input and representation releases bound by a comparison receipt."""

    def __init__(self, receipt, inputs, *, representation=None, representation_sha256=None):
        self.data, self.table, self.provenance = _bind(inputs, representation, representation_sha256)
        if json.loads(Path(receipt).read_text()) != self.provenance:
            raise ValueError("representation comparison receipt differs from verified inputs")

    def features(self, arm, rows):
        parents = self.data.obs.iloc[rows].drug_parent_inchikey
        return self.table.vectors(arm, parents)


def evaluate(comparison, arm, fold, predictions, prediction_ids, output):
    """Score a complete cold-drug fold and retain its parent-equal study statistic."""
    if arm not in (*ARMS, "Mean_cell") or fold not in range(3):
        raise ValueError("unknown arm or cold-drug fold")
    data = comparison.data
    samples, sample_metrics = score(data, SCENARIO, fold, "test", predictions, prediction_ids)
    parent_scores = samples.groupby("drug_parent", as_index=False).pcc.mean()
    provenance = comparison.provenance
    generation = provenance["view_generation"]
    summary = {"source": data.source, "fold": fold, "arm": arm, "scenario": SCENARIO,
               "seed": provenance["seed"], "pair_sha256": provenance["pair_sha256"],
               "view_sha256": generation["view_sha256"], "axis_sha256": generation["axis_sha256"],
               "corpus_sha256": generation["corpus_sha256"],
               "representation_sha256": provenance["representation_sha256"],
               "parent_equal_pcc": parent_scores.pcc.mean(),
               "defined_parents": int(parent_scores.pcc.count()),
               "total_parents": len(parent_scores), **sample_metrics}

    output = Path(output)
    if output.exists():
        raise FileExistsError(f"evaluation output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".evaluate-", dir=output.parent) as temporary:
        candidate = Path(temporary) / "report"
        candidate.mkdir()
        samples.to_csv(candidate / "per_sample.csv", index=False)
        parent_scores.to_csv(candidate / "parent_pcc.csv", index=False)
        pd.DataFrame([summary]).to_csv(candidate / "summary.csv", index=False)
        candidate.rename(output)
    return summary


def render(summary_csv, output):
    """Render the parent-equal Delta PCC comparison from recorded CSV values."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from stage2_util import palette, plotting

    rows = plotting.read_rows(summary_csv)
    if not rows:
        raise ValueError("representation summary CSV is empty")
    figure = Figure(figsize=(8, 4), layout="constrained")
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    colors = {"tenv": "blue", "tahoe_bulk": "orange"}
    for index, row in enumerate(rows):
        family = "gray" if row["arm"] == "Mean_cell" else colors[row["source"]]
        value = float(row["parent_equal_pcc"]) if row["parent_equal_pcc"] else np.nan
        if np.isfinite(value):
            axis.bar(index, value, color=palette.color(family, palette.FILL_SHADE),
                     edgecolor=palette.color(family, palette.LINE_SHADE))
        else:
            axis.scatter(index, 0, marker="x", color=palette.color("gray", palette.REFERENCE_SHADE))
    labels = [f"{row['arm']}\n{row['source']} · fold {row['fold']}" for row in rows]
    axis.set_xticks(range(len(rows)), labels, rotation=25, ha="right")
    axis.set_ylabel("Parent-equal Δ Pearson")
    axis.set_axisbelow(True)
    axis.grid(axis="y", color=palette.color("gray", 8))
    plotting.save_figure(figure, output, dpi=180)
