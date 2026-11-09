"""Bounded reconstruction fixtures exercise the fit and replay workflow."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "script"
sys.path.insert(0, str(SCRIPT))


class BenchmarkTest(unittest.TestCase):
    def setUp(self):
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.root = Path(work.name)
        self.release = self.root / "release"
        self.release.mkdir()
        self.split_path = self.root / "split.json"
        self.artifact = self.root / "mean"
        self.profile_plates = {"tahoe": ["A", "B", "C", "A"], "tenv": ["A", "B", "D"]}
        self.expressions = {
            "tahoe": np.array([[1, 2, 3], [2, 4, 8], [3, 1, 2], [4, 6, 8]], dtype=np.float32),
            "tenv": np.array([[8, 3, 0], [5, 5, 5], [0, 6, 4]], dtype=np.float32),
        }
        (self.release / "genes.csv").write_text("symbol,hgnc_id\nX,HGNC:1\nY,HGNC:2\nZ,HGNC:3\n")
        for source, expression in self.expressions.items():
            np.save(self.release / f"{source}_expression.npy", expression)
            with (self.release / f"{source}_profiles.csv").open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["source_row", "profile_id", "plate", "is_control"])
                for index, plate in enumerate(self.profile_plates[source]):
                    writer.writerow([index * 2, f"{source}-{index}", plate, index == 0])
            (self.release / f"{source}_mapping.csv").write_text("fixture mapping\n")
        self.write_manifest()

    def write_manifest(self):
        files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in self.release.iterdir() if p.name != "manifest.json"}
        manifest = {
            "schema_version": 1,
            "fixture": True,
            "normalization": "log1p(CP10K of raw counts projected to the ordered CM axis)",
            "inputs": {"genes.csv": files["genes.csv"]},
            "sources": {source: {"shape": list(expr.shape)}
                        for source, expr in self.expressions.items()},
            "files": files,
        }
        (self.release / "manifest.json").write_text(json.dumps(manifest))

    def cli(self, action, *args, success=True, script="reconstruct.py"):
        result = subprocess.run(
            [sys.executable, str(SCRIPT / script), action, "--fixture",
             "--release", str(self.release), "--split", str(self.split_path), *map(str, args)],
            capture_output=True, text=True,
        )
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
            return json.loads(result.stdout)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        return result.stderr

    def test_mean_fit_and_independent_reload_use_frozen_roles(self):
        from benchmark import Corpus, Split

        split_record = self.cli("freeze-split")
        split = Split(Corpus(self.release, fixture=True), self.split_path)
        training = []
        for source, plates in self.profile_plates.items():
            rank = sorted(set(plates), key=lambda plate: hashlib.sha256(
                json.dumps([20261109, source, plate]).encode()).digest())
            expected = {"validation": rank[:1], "test": rank[1:2], "train": rank[2:]}
            for role, role_plates in expected.items():
                rows = [i for i, plate in enumerate(plates) if plate in role_plates]
                np.testing.assert_array_equal(split.rows[source][role], rows)
                record = split_record["sources"][source]["roles"][role]
                self.assertEqual(record["profiles"], len(rows))
                self.assertEqual(record["plates"], role_plates)
                self.assertEqual(len(record["membership_sha256"]), 64)
                if role == "train":
                    training.extend(self.expressions[source][rows])
        self.cli("fit", "--output", self.artifact, "--batch-size", 1)
        mean = np.load(self.artifact / "mean.npy")
        np.testing.assert_array_equal(mean, np.mean(np.array(training), axis=0))

        for role in ("train", "validation", "test"):
            path = self.root / f"{role}.json"
            args = ["--artifact", self.artifact, "--role", role, "--batch-size", 1]
            scores = self.cli("score", *args, "--output", path)
            replay = self.cli("reload", *args, "--expected", path)
            self.assertEqual(scores, replay)
            for source, result in scores["sources"].items():
                self.assertEqual(result["profiles"], len(split.rows[source][role]))
                self.assertEqual(result["relative_squared_error"], 1.0)

        with (self.artifact / "mean.npy").open("ab") as handle:
            handle.write(b"changed")
        self.assertIn("digest", self.cli("reload", "--artifact", self.artifact, "--expected", path, success=False))

    def test_identity_and_undefined_metrics(self):
        from benchmark import Corpus, Split, score

        self.cli("freeze-split")
        corpus = Corpus(self.release, fixture=True)
        split = Split(corpus, self.split_path)
        for role in ("train", "validation", "test"):
            scores = score(corpus, split, role, lambda values: values, np.zeros(3), 1)
            for source, result in scores["sources"].items():
                rows = split.rows[source][role]
                defined = int(np.count_nonzero(np.ptp(self.expressions[source][rows], axis=1)))
                self.assertEqual(result["profiles"], len(rows))
                self.assertEqual(result["defined_pcc_profiles"], defined)
                if defined:
                    self.assertAlmostEqual(result["pcc"], 1.0)
                else:
                    self.assertIsNone(result["pcc"])
                self.assertEqual(result["rmse"], 0.0)
                self.assertEqual(result["relative_squared_error"], 0.0)
        role = next(role for role, rows in split.rows["tenv"].items() if 1 in rows)
        scores = score(corpus, split, role, lambda values: values, np.full(3, 5), 2)
        self.assertIsNone(scores["sources"]["tenv"]["relative_squared_error"])

    def prepare_pca_fixture(self):
        genes = 577
        center = np.linspace(2, 3, genes, dtype=np.float32)
        offsets = np.zeros(genes, dtype=np.float32)
        offsets[[0, 31, 32, 127, 128, 511, 512]] = [1, 2, 3, 4, 5, 6, 7]
        (self.release / "genes.csv").write_text(
            "symbol,hgnc_id\n" + "".join(f"G{i},HGNC:{i + 1}\n" for i in range(genes))
        )
        # Diagonal training covariance with eigenvalue gaps at the three cutoffs.
        axes = np.diag(np.repeat(np.array([1, 0.75, 0.5, 0.25], dtype=np.float32),
                                 [32, 96, 384, genes - 512]))
        half = len(axes) // 2
        for source, scale, directions in (("tahoe", 1, axes[:half]), ("tenv", 2, axes[half:])):
            ranked = sorted("ABC", key=lambda plate: hashlib.sha256(
                json.dumps([20261109, source, plate]).encode()).digest())
            values = np.vstack((center + scale * offsets, center - offsets,
                                center + directions, center - directions))
            self.expressions[source] = values
            np.save(self.release / f"{source}_expression.npy", values)
            with (self.release / f"{source}_profiles.csv").open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["source_row", "profile_id", "plate", "is_control"])
                for index in range(len(values)):
                    writer.writerow([index, f"{source}-{index}", ranked[min(index, 2)], False])
        self.write_manifest()
        return center, offsets

    def test_pca_fit_and_independent_reload_at_each_rank(self):
        center, offsets = self.prepare_pca_fixture()
        genes = len(center)
        self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        models = self.root / "pca"
        self.cli("fit", "--model", "PCA", "--output", models)
        for rank in (32, 128, 512):
            with self.subTest(rank=rank):
                model_id = f"PCA_{rank}"
                artifact = models / model_id
                self.assertEqual(np.load(artifact / "components.npy").shape, (rank, genes))
                np.testing.assert_allclose(np.load(artifact / "mean.npy"), center, atol=1e-7)
                args = ["--artifact", artifact, "--mean-artifact", self.artifact]
                path = self.root / f"{model_id}.json"
                scores = self.cli("score", *args, "--output", path)
                self.assertEqual(scores["model_id"], model_id)
                self.assertEqual(scores, self.cli("reload", *args, "--expected", path))
                for source, scale in (("tahoe", 1), ("tenv", 2)):
                    result = scores["sources"][source]
                    self.assertEqual(result["profiles"], 1)
                    expected_mse = float(np.square(offsets[rank:]).sum()) * scale ** 2 / genes
                    self.assertAlmostEqual(result["mse"], expected_mse, places=6)

    def test_comparison_seals_validation_and_replays_all_test_models(self):
        self.prepare_pca_fixture()
        self.cli("freeze-split")
        models = self.root / "models"
        seal_path = self.root / "seal.json"
        score_path = self.root / "scores.json"
        args = ("--models", models, "--seal", seal_path)
        self.assertIn(str(seal_path), self.cli("score", *args, "--output", score_path,
                                              success=False, script="compare.py"))
        self.cli("fit", "--model", "PCA", "--output", models)
        model_ids = ("mean", "PCA_32", "PCA_128", "PCA_512", "MLP_32", "MLP_128", "MLP_512", "AAE_512")
        for model in model_ids:
            artifact = models / model
            if not model.startswith("PCA_"):
                self.cli("fit", "--model", model, "--output", artifact)
            self.cli("score", "--artifact", artifact, "--mean-artifact", models / "mean",
                     "--output", self.root / f"{model}-validation.json")

        validation_path = self.root / "AAE_512-validation.json"
        validation = validation_path.read_text()
        validation_path.write_text(validation.replace('"role": "validation"', '"role": "test"'))
        self.assertIn("validation", self.cli("seal", *args, "--validation-dir", self.root,
                                             success=False, script="compare.py"))
        self.assertFalse(seal_path.exists())
        validation_path.write_text(validation)
        seal = self.cli("seal", *args, "--validation-dir", self.root, script="compare.py")
        self.assertEqual(set(seal["models"]), set(model_ids))
        scores = self.cli("score", *args, "--output", score_path, script="compare.py")
        self.assertEqual(scores["role"], "test")
        self.assertEqual(scores["seal_sha256"], hashlib.sha256(seal_path.read_bytes()).hexdigest())
        self.assertEqual(scores["models"]["AAE_512"]["latent"], {"tokens": 512, "scalars": 65536})
        self.assertEqual(scores["models"]["PCA_512"]["latent"], {"tokens": None, "scalars": 512})
        self.assertEqual(set(scores["models"]), set(model_ids))
        for model in model_ids:
            for source in ("tahoe", "tenv"):
                self.assertEqual(scores["models"][model]["sources"][source]["profiles"], 1)
        self.assertEqual(scores, self.cli("reload", *args, "--expected", score_path, script="compare.py"))

        figure_path = self.root / "comparison.png"
        figure = subprocess.run(
            [sys.executable, str(SCRIPT / "compare.py"), "figure", "--expected", str(score_path),
             "--output", str(figure_path)], capture_output=True, text=True,
        )
        self.assertEqual(figure.returncode, 0, figure.stderr)
        self.assertTrue(figure_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

        manifest_path = models / "mean" / "manifest.json"
        manifest_path.write_text(manifest_path.read_text() + "\n")
        self.assertIn("seal", self.cli("reload", *args, "--expected", score_path,
                                       success=False, script="compare.py"))

    def test_mlp_reloads_the_tenv_selected_checkpoint_at_each_latent_width(self):
        for source, values in self.expressions.items():
            validation_plate = min(self.profile_plates[source], key=lambda plate: hashlib.sha256(
                json.dumps([20261109, source, plate]).encode()).digest())
            for index, plate in enumerate(self.profile_plates[source]):
                # Training moves outputs towards four and away from tenv validation.
                values[index] = {"tenv": 0, "tahoe": 6}[source] if plate == validation_plate else 4
            np.save(self.release / f"{source}_expression.npy", values)
        self.write_manifest()
        split = self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        for width in (32, 128, 512):
            with self.subTest(width=width):
                model_id = f"MLP_{width}"
                artifact = self.root / model_id
                result = self.cli("fit", "--model", model_id, "--output", artifact)
                metadata = result["metadata"]
                self.assertEqual(metadata["latent_scalars"], width)
                self.assertIsNone(metadata["latent_tokens"])
                self.assertEqual(metadata["fit_profiles"], {
                    source: split["sources"][source]["roles"]["train"]["profiles"]
                    for source in self.expressions
                })
                history = json.loads((artifact / "epochs.json").read_text())
                self.assertEqual(metadata["selected_epoch"], 1)
                self.assertLess(history[0]["validation_mse"]["tenv"],
                                history[-1]["validation_mse"]["tenv"])
                if width == 32:
                    self.assertGreater(history[0]["validation_mse"]["tahoe"],
                                       history[-1]["validation_mse"]["tahoe"])
                expected = self.root / f"{model_id}.json"
                args = ["--artifact", artifact, "--mean-artifact", self.artifact,
                        "--role", "validation", "--batch-size", 2]
                score = self.cli("score", *args, "--output", expected)
                self.assertEqual(score["provenance"], result["provenance"])
                for source in self.expressions:
                    self.assertAlmostEqual(
                        score["sources"][source]["mse"], history[0]["validation_mse"][source],
                    )
                self.assertEqual(score, self.cli("reload", *args, "--expected", expected))

    def test_aae_replays_the_tenv_selected_full_axis_checkpoint(self):
        split = self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        model_id = "AAE_512"
        artifact = self.root / model_id
        result = self.cli("fit", "--model", model_id, "--output", artifact)
        metadata = result["metadata"]
        self.assertEqual(result["model_id"], model_id)
        self.assertEqual(metadata["latent_tokens"], 512)
        self.assertEqual(metadata["latent_scalars"], 65_536)
        self.assertEqual(metadata["architecture"]["width"], 128)
        self.assertEqual(metadata["architecture"]["heads"], 8)
        self.assertEqual(metadata["fit_profiles"], {
            source: split["sources"][source]["roles"]["train"]["profiles"]
            for source in self.expressions
        })
        history = json.loads((artifact / "epochs.json").read_text())
        selected = min(history, key=lambda epoch: (epoch["validation_mse"]["tenv"], epoch["epoch"]))
        self.assertEqual(metadata["selected_epoch"], selected["epoch"])
        expected = self.root / "aae-validation.json"
        args = ["--artifact", artifact, "--mean-artifact", self.artifact, "--role", "validation"]
        scores = self.cli("score", *args, "--output", expected)
        self.assertEqual(scores["provenance"], result["provenance"])
        for source, score in scores["sources"].items():
            self.assertEqual(score["profiles"], split["sources"][source]["roles"]["validation"]["profiles"])
            self.assertAlmostEqual(score["mse"], selected["validation_mse"][source])
        self.assertEqual(scores, self.cli("reload", *args, "--expected", expected))
        self.assertIn("frozen batch size", self.cli("score", *args, "--batch-size", 32, success=False))

    def test_refuses_changed_split_and_release(self):
        self.cli("freeze-split")
        frozen = json.loads(self.split_path.read_text())
        frozen["sources"]["tahoe"]["roles"]["test"]["plates"] = ["other"]
        self.split_path.write_text(json.dumps(frozen))
        self.assertIn("split", self.cli("fit", "--output", self.artifact, success=False))
        with (self.release / "tenv_expression.npy").open("ab") as handle:
            handle.write(b"changed")
        self.assertIn("digest", self.cli("fit", "--output", self.artifact, success=False))


if __name__ == "__main__":
    unittest.main()
