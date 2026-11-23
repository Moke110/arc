import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from anndata import AnnData
import h5py
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "output/api"))
from tenv5k_raw_first import workflow
from expression_view import PAYLOAD_SCHEMA, RECIPE, SCHEMA
from prediction_axis import open_axes
from stage2_util import metrics, palette, plotting
from stage2_util.representations import ARMS, SCHEMA as REPRESENTATIONS, WIDTHS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixtures(root):
    view = root / "view"
    view.mkdir(parents=True)
    axis = open_axes().axis(5000)
    obs = pd.DataFrame(
        [("vehicle", "", "", 0), *[("treated", f"p{parent}", f"smiles{parent}", 1)
                                     for parent in range(6) for _ in range(3)],
         ("treated", "p0", "smiles0", 10), ("vehicle", "", "", 0)],
        columns=["role", "drug_parent_inchikey", "drug_parent_smiles", "dose_um"],
        index=[f"well-{i}" for i in range(21)]).assign(cell_line="hiPSC-CM", matched_control_group="plate")
    obs["profile_id"] = obs.index
    var = pd.DataFrame({"gene_index": axis.gene_indices, "symbol": axis.symbols},
                       index=list(map(str, axis.gene_indices)))
    expression = np.zeros((len(obs), 5000), dtype=np.float32)
    expression[:, :3] = [5, 10, 20]
    expression[0, :3] = [8, 9, 10]
    expression[-1, :3] = [4, 5, 6]
    control = np.zeros_like(expression)
    control[:, :3] = [1, 2, 3]
    source = view / "tenv.h5ad"
    AnnData(expression, obs=obs, var=var,
            layers={"control": control, "delta": expression - control}).write_h5ad(source)
    tags = {"schema": SCHEMA, "axis_sha256": axis.sha256, "corpus_sha256": "a" * 64}
    with h5py.File(source, "r+") as handle:
        handle.attrs.update(dict(tags, schema=PAYLOAD_SCHEMA), source="tenv", family="full")
    exclusions = view / "tenv.exclusions.csv"
    exclusions.write_text("profile_id,role,reason\n")
    membership = view / "tenv.members.csv"
    membership.write_text("profile_id,member_profile_id,role\n")
    (view / "manifest.json").write_text(json.dumps(dict(tags, **RECIPE, axis_size=5000, sources=["tenv"],
        outputs={str(p.relative_to(view)): digest(p) for p in (source, exclusions, membership)})))

    vectors = root / "representations"
    vectors.mkdir()
    parents = pd.DataFrame({"drug_parent_inchikey": [f"p{i}" for i in reversed(range(6))],
                            "drug_parent_smiles": [f"smiles{i}" for i in reversed(range(6))]})
    parents.to_csv(vectors / "parents.csv", index=False)
    parents[["drug_parent_inchikey"]].assign(conformers=1).to_csv(vectors / "e3fp_conformers.csv", index=False)
    for arm in ARMS:
        values = np.zeros((6, WIDTHS[arm]), dtype=np.float32)
        values[:, 0] = [1, 0] * 3
        np.save(vectors / f"{arm}.npy", values)
    manifest = vectors / "manifest.json"
    manifest.write_text(json.dumps({"schema": REPRESENTATIONS, "source": "tenv", "widths": WIDTHS,
                                    "inputs_sha256": "c" * 64,
                                    "files": {p.name: digest(p) for p in vectors.iterdir()}}))
    return view, vectors, digest(manifest)


class RawFirstSingleDoseTest(unittest.TestCase):
    def test_preserves_controls_and_shared_inputs_then_scores_and_renders(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            view, vectors, expected = fixtures(root)
            directory = root / "inputs"
            workflow.prepare(directory, view_path=view, representations=vectors,
                             representation_sha256=expected)
            data = workflow.Inputs(directory)
            self.assertEqual(list(data.obs.index), [f"well-{i}" for i in range(1, 19)])
            self.assertEqual(data.controls[0], [18, 19])
            self.assertEqual(data.expression[18:20, :3].tolist(), [[8, 9, 10], [4, 5, 6]])
            self.assertEqual(data.control[0, :3].tolist(), [1, 2, 3])
            self.assertEqual(data.delta[0, :3].tolist(), [4, 8, 17])
            self.assertEqual(data.fingerprints([0, 3])[:, 0].tolist(), [0, 1])
            self.assertEqual(data.config["training_seeds"], [17, 42, 2026])
            rows = data.rows("warm-start", 0, "test")
            prediction = root / "prediction.npy"
            np.save(prediction, data.expression[rows])
            identities = root / "ids.csv"
            pd.DataFrame({"profile_id": data.obs.index[rows]}).to_csv(identities, index=False)
            report = root / "report"
            with patch.object(metrics, "vector_metrics", wraps=metrics.vector_metrics) as calls:
                workflow.evaluate(directory, prediction, identities, report, scenario="warm-start", fold=0, model="Mean", seed=42, run_id="mean-42")
                self.assertGreater(calls.call_count, 2 * len(rows))
            summary = pd.read_csv(report / "summary.csv")
            self.assertEqual(summary.loc[0, ["model", "seed", "run_id"]].tolist(), ["Mean", 42, "mean-42"])
            per_parent = pd.read_csv(report / "per_parent.csv")
            self.assertTrue(per_parent.rmse.eq(0).all())
            self.assertEqual(per_parent.n_wells.sum(), len(rows))
            comparison = root / "comparison.csv"
            pd.concat([summary, summary.assign(model="PRnet", run_id="prnet-42")]).to_csv(comparison, index=False)
            image = root / "figure.png"
            with (patch.object(plotting, "read_rows", wraps=plotting.read_rows) as read,
                  patch.object(palette, "color", wraps=palette.color) as color,
                  patch.object(plotting, "save_figure", wraps=plotting.save_figure) as save):
                workflow.render(comparison, image)
                read.assert_called_once()
                self.assertGreater(color.call_count, 0)
                figure = save.call_args.args[0]
                self.assertEqual(len(figure.axes[0].patches), 2)
                labels = [tick.get_text() for tick in figure.axes[0].get_xticklabels()]
                self.assertIn("Mean · 42 · mean-42", labels[0])
                self.assertIn("PRnet · 42 · prnet-42", labels[1])
            self.assertTrue(image.read_bytes().startswith(b"\x89PNG"))
