import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

from anndata import AnnData
import h5py
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "output/api"))
from tsg_raw_first import preparation
from tsg_raw_first.inputs import Inputs


class RawFirstInputsTest(unittest.TestCase):
    def test_pair_keeps_full_view_control_and_individual_vehicles(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            views = make_views(root / "view")
            representations = make_representations(root)
            pair = root / "inputs"
            preparation.prepare(pair, view=views, representations=representations)

            for source in ("tenv", "tahoe_bulk"):
                with self.subTest(source=source):
                    data = Inputs(pair / source)
                    self.assertEqual(data.obs.index[0], f"{source}-treated-0-0")
                    self.assertEqual(data.controls[0], [18, 19])
                    self.assertEqual(data.expression[0, :3].tolist(), [2., 4., 8.])
                    self.assertEqual(data.control[0, :3].tolist(), [1., 2., 3.])
                    self.assertEqual(data.expression[18:20, :3].tolist(), [[5., 6., 7.], [1., 1., 2.]])
                    self.assertEqual(data.fingerprints([0, 3])[:, :2].tolist(), [[1., 0.], [0., 1.]])
                    self.assertEqual(data.config["view_generation"]["view_sha256"],
                                     sha256(views / "manifest.json"))

            source_manifest = pair / "tahoe_bulk/manifest.json"
            source_manifest.write_bytes((pair / "tenv/manifest.json").read_bytes())
            with self.assertRaisesRegex(ValueError, "outside the pinned pair"):
                Inputs(pair / "tahoe_bulk")

    def test_kpgt_wrapper_aligns_vectors_to_the_new_pair(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            representations = make_representations(root)
            pair = root / "inputs"
            preparation.prepare(pair, view=make_views(root / "view"), representations=representations)
            destination = root / "kpgt"
            from tsg_raw_first.inputs import KPGTInputs
            preparation.prepare_kpgt(pair / "tenv", destination, representation=representations["tenv"])
            data = KPGTInputs(destination)
            self.assertEqual(data.fingerprints([0, 3])[:, 0].tolist(), [-4., -3.])
            self.assertEqual(data.pair_sha256, sha256(pair / "manifest.json"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_manifest(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def make_views(root):
    from prediction_axis import open_axes
    from expression_view import PAYLOAD_SCHEMA, RECIPE, SCHEMA, SOURCES

    axis = open_axes().axis(5000)
    var = pd.DataFrame({"gene_index": axis.gene_indices, "symbol": axis.symbols},
                       index=[str(gene) for gene in axis.gene_indices])
    root.mkdir(parents=True)
    manifest = {"schema": SCHEMA, "axis_sha256": axis.sha256, "axis_size": 5000,
                "corpus_sha256": "a" * 64, "sources": list(SOURCES), **RECIPE, "outputs": {}}
    for source in SOURCES:
        ids = [f"{source}-vehicle-0"]
        rows = [("vehicle", "", "", 0.)]
        doses = (0.2, 1., 10.) if source == "tenv" else (0.05, 0.5, 5.)
        for parent in range(6):
            for dose, value in enumerate(doses):
                ids.append(f"{source}-treated-{parent}-{dose}")
                rows.append(("treated", f"parent-{parent}", f"SMILES-{parent}", value))
        ids.append(f"{source}-vehicle-1")
        rows.append(("vehicle", "", "", 0.))
        obs = pd.DataFrame(rows, columns=["role", "drug_parent_inchikey", "drug_parent_smiles", "dose_um"],
                           index=ids).assign(profile_id=ids, cell_line="line", matched_control_group="plate")
        matrix = np.zeros((len(obs), 5000), dtype=np.float32)
        matrix[:, :3] = [2, 4, 8]
        matrix[0, :3] = [5, 6, 7]
        matrix[-1, :3] = [1, 1, 2]
        control = np.zeros_like(matrix)
        control[:, :3] = [1, 2, 3]
        payload = f"{source}.h5ad"
        AnnData(matrix, obs=obs, var=var, layers={"control": control, "delta": matrix-control}).write_h5ad(
            root / payload
        )
        with h5py.File(root / payload, "r+") as handle:
            handle.attrs["schema"] = PAYLOAD_SCHEMA
            for key in ("axis_sha256", "corpus_sha256"):
                handle.attrs[key] = manifest[key]
            handle.attrs["source"] = source
            handle.attrs["family"] = "full"
        exclusion = f"{source}.exclusions.csv"
        (root / exclusion).write_text("profile_id,role,reason\n")
        membership = f"{source}.members.csv"
        (root / membership).write_text("profile_id,member_profile_id,role\n")
        manifest["outputs"].update({name: sha256(root / name) for name in (payload, exclusion, membership)})
    write_manifest(root / "manifest.json", manifest)
    return root


def make_representations(root):
    from stage2_util.representations import ARMS, SCHEMA, WIDTHS

    result = {}
    for source in ("tenv", "tahoe_bulk"):
        out = root / f"representations-{source}"
        out.mkdir()
        order = list(reversed(range(6)))
        parents = pd.DataFrame({"drug_parent_inchikey": [f"parent-{i}" for i in order],
                                "drug_parent_smiles": [f"SMILES-{i}" for i in order]})
        parents.to_csv(out / "parents.csv", index=False)
        parents[["drug_parent_inchikey"]].assign(conformers=1).to_csv(
            out / "e3fp_conformers.csv", index=False
        )
        for arm in ARMS:
            data = np.zeros((6, WIDTHS[arm]), dtype=np.float32)
            if arm == "KPGT":
                data[:, 0] = np.asarray(order) - 4
            else:
                data[np.arange(6), np.asarray(order) % 2] = 1
            np.save(out / f"{arm}.npy", data)
        outputs = {p.name: sha256(p) for p in out.iterdir()}
        write_manifest(out / "manifest.json", {"schema": SCHEMA, "source": source,
                       "widths": WIDTHS, "inputs_sha256": "b" * 64, "files": outputs})
        result[source] = (out, sha256(out / "manifest.json"))
    return result


class RawFirstEvaluationTest(unittest.TestCase):
    def test_evaluate_render_and_reject_reordered_prediction_ids(self):
        from unittest.mock import patch
        sys.path.insert(0, str(STEP / "execution"))
        import raw_first
        from stage2_util import metrics, palette, plotting

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pair = root / "inputs"
            preparation.prepare(pair, view=make_views(root / "view"), representations=make_representations(root))
            data = Inputs(pair / "tenv")
            rows = data.rows("cold-drug", 0, "test")
            expected = data.expression[rows].copy()
            prediction = root / "prediction.npy"
            np.save(prediction, expected)
            ids = root / "ids.csv"
            pd.DataFrame({"profile_id": data.obs.index[rows]}).to_csv(ids, index=False)
            report = root / "report"
            args = ["evaluate", "--inputs", str(data.directory), "--scenario", "cold-drug", "--fold", "0",
                    "--method", "no-dose", "--seed", "42", "--run-id", "no-dose-42",
                    "--prediction", str(prediction), "--prediction-ids", str(ids), "--output", str(report)]
            with patch.object(metrics, "vector_metrics", wraps=metrics.vector_metrics) as called:
                raw_first.main(args)
                self.assertEqual(called.call_count, len(rows))
            first_truth, first_predicted = called.call_args_list[0].args
            self.assertEqual(first_truth[:3].tolist(), [1., 2., 5.])
            self.assertTrue(np.array_equal(first_truth, first_predicted))
            summary = pd.read_csv(report / "summary.csv")
            self.assertEqual(summary.loc[0, ["method", "seed", "run_id"]].tolist(), ["no-dose", 42, "no-dose-42"])
            self.assertEqual(summary.iloc[0].rmse, 0.)
            self.assertEqual(summary.iloc[0].r2, 1.)
            self.assertEqual(summary.iloc[0].n_samples, len(rows))

            csv_only = root / "minimal.csv"
            csv_only.write_bytes((report / "summary.csv").read_bytes())
            png = root / "figures/metrics.png"
            with (patch.object(plotting, "read_rows", wraps=plotting.read_rows) as read,
                  patch.object(plotting, "save_figure", wraps=plotting.save_figure) as save,
                  patch.object(palette, "color", wraps=palette.color) as color):
                raw_first.main(["render", "--summary", str(csv_only), "--output", str(png)])
                read.assert_called_once_with(csv_only)
                save.assert_called_once()
                self.assertGreater(color.call_count, 1)
            self.assertTrue(png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

            pd.DataFrame({"profile_id": data.obs.index[rows][::-1]}).to_csv(ids, index=False)
            args[-1] = str(root / "bad-report")
            with self.assertRaisesRegex(ValueError, "prediction identity"):
                raw_first.main(args)
            self.assertFalse((root / "bad-report").exists())
