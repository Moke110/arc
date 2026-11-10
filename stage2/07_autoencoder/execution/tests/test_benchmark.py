"""Bounded reconstruction fixtures exercise the fit and replay workflow."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

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

    def assert_comparison_png(self, action, input_flag, input_path, output_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT / "compare.py"), action, input_flag, str(input_path),
             "--output", str(output_path)], capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(output_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

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

    def test_only_comparison_cli_accepts_real_test_scoring(self):
        for action, arguments in (("score", []), ("reload", ["--expected", "unused.json"])):
            with self.subTest(action=action):
                result = subprocess.run(
                    [sys.executable, str(SCRIPT / "reconstruct.py"), action, "--role", "test",
                     "--allow-test", "--release", str(self.root / "absent"), *arguments],
                    capture_output=True, text=True,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("real test scoring requires a selection seal through compare.py",
                              result.stderr)

    def test_batch_limit_applies_to_cli_and_corpus_reads(self):
        from benchmark import Corpus, Split, batches

        self.cli("freeze-split")
        self.assertIn("between 1 and 1024", self.cli("fit", "--output", self.artifact,
                                                   "--batch-size", 1025, success=False))
        corpus = Corpus(self.release, fixture=True)
        split = Split(corpus, self.split_path)
        with self.assertRaisesRegex(ValueError, "between 1 and 1024"):
            list(batches(corpus, split, "tahoe", "train", 1025))
        np.testing.assert_array_equal(corpus.read("tahoe", np.zeros(1024, dtype=int)),
                                      np.tile(self.expressions["tahoe"][0], (1024, 1)))
        with self.assertRaisesRegex(ValueError, "at most 1024 rows"):
            corpus.read("tahoe", np.zeros(1025, dtype=int))

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

    def test_expanded_comparison_requires_a_separate_seal_path(self):
        error = self.cli("seal", "--expanded-aae", success=False, script="compare.py")
        self.assertIn("expanded comparison requires a separate --seal path", error)

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
        model_ids = ("mean", "PCA_32", "PCA_128", "PCA_512", "MLP_32", "MLP_128", "MLP_512", "AAE_512",
                     "MatrixAE_64x8")
        for model in model_ids:
            artifact = models / model
            if not model.startswith("PCA_"):
                self.cli("fit", "--model", model, "--output", artifact)
            self.cli("score", "--artifact", artifact, "--mean-artifact", models / "mean",
                     "--output", self.root / f"{model}-validation.json")

        validation_path = self.root / "MatrixAE_64x8-validation.json"
        validation = validation_path.read_text()
        validation_path.write_text(validation.replace('"role": "validation"', '"role": "test"'))
        self.assertIn("validation", self.cli("seal", *args, "--validation-dir", self.root,
                                             success=False, script="compare.py"))
        self.assertFalse(seal_path.exists())
        validation_path.write_text(validation)
        seal = self.cli("seal", *args, "--validation-dir", self.root, script="compare.py")
        self.assertEqual(set(seal["models"]), set(model_ids))
        matrix = seal["models"]["MatrixAE_64x8"]
        history = json.loads((models / "MatrixAE_64x8" / "epochs.json").read_text())
        best = min(history, key=lambda check: (check["validation_mse"]["tenv"], check["epoch"]))
        self.assertEqual(matrix["selected_epoch"], best["epoch"])
        self.assertEqual(matrix["device"], "cpu")
        scores = self.cli("score", *args, "--output", score_path, script="compare.py")
        self.assertEqual(scores["role"], "test")
        self.assertEqual(scores["seal_sha256"], hashlib.sha256(seal_path.read_bytes()).hexdigest())
        self.assertEqual(scores["models"]["AAE_512"]["latent"], {"tokens": 512, "scalars": 65536})
        self.assertEqual(scores["models"]["PCA_512"]["latent"], {"tokens": None, "scalars": 512})
        self.assertEqual(scores["models"]["MatrixAE_64x8"]["latent"], {"tokens": 64, "scalars": 512})
        self.assertEqual(set(scores["models"]), set(model_ids))
        for model in model_ids:
            for source in ("tahoe", "tenv"):
                self.assertEqual(scores["models"][model]["sources"][source]["profiles"], 1)
        self.assertEqual(scores, self.cli("reload", *args, "--expected", score_path, script="compare.py"))

        self.assert_comparison_png("figure", "--expected", score_path, self.root / "comparison.png")

        legacy_seal_path = self.root / "legacy-seal.json"
        legacy_score_path = self.root / "legacy-scores.json"
        del seal["models"]["MatrixAE_64x8"]
        legacy_seal_path.write_text(json.dumps(seal))
        legacy_args = ("--models", models, "--seal", legacy_seal_path)
        legacy_scores = self.cli("score", *legacy_args, "--output", legacy_score_path, script="compare.py")
        self.assertEqual(legacy_scores["models"], {
            model: scores["models"][model] for model in model_ids[:-1]
        })
        self.assertEqual(
            legacy_scores,
            self.cli("reload", *legacy_args, "--expected", legacy_score_path, script="compare.py"),
        )
        self.assert_comparison_png(
            "figure", "--expected", legacy_score_path, self.root / "legacy-comparison.png",
        )

        self.assertEqual(seal["schema_version"], 1)
        self.assertEqual(scores["schema_version"], 1)
        self.assertNotIn("display_name", scores["models"]["AAE_512"])
        self.assert_expanded_comparison(models, scores)

        manifest_path = models / "mean" / "manifest.json"
        manifest_path.write_text(manifest_path.read_text() + "\n")
        self.assertIn("seal", self.cli("reload", *args, "--expected", score_path,
                                       success=False, script="compare.py"))
        self.assertIn("seal", self.cli("reload", *legacy_args, "--expected", legacy_score_path,
                                       success=False, script="compare.py"))

        del seal["models"]["AAE_512"]
        legacy_seal_path.write_text(json.dumps(seal))
        self.assertIn("unexpected comparison model set", self.cli(
            "reload", *legacy_args, "--expected", legacy_score_path, success=False, script="compare.py",
        ))

    def assert_expanded_comparison(self, models, historical_scores):
        seal_path = self.root / "expanded-seal.json"
        score_path = self.root / "expanded-scores.json"
        args = ("--models", models, "--seal", seal_path)
        seal_args = (*args, "--expanded-aae", "--validation-dir", self.root)
        variants = ("AAE_512_Transformer", "AAE_512+16")
        for model in variants:
            artifact = models / model
            self.cli("fit", "--model", model, "--output", artifact)
            self.cli("score", "--artifact", artifact, "--mean-artifact", models / "mean",
                     "--output", self.root / f"{model}-validation.json")

        self.assertIn(str(seal_path), self.cli(
            "score", *args, "--output", score_path, success=False, script="compare.py",
        ))
        self.assertFalse(score_path.exists())
        self.assert_invalid_expanded_selection(models, seal_path, seal_args)
        seal = self.cli("seal", *seal_args, script="compare.py")
        scores = self.cli("score", *args, "--output", score_path, script="compare.py")
        self.assertEqual(seal["schema_version"], 2)
        self.assertEqual(scores["schema_version"], 2)
        self.assertEqual(set(scores["models"]), {*historical_scores["models"], *variants})
        for artifact_id, display_name, tokens, scalars in (
            ("AAE_512", "AAE_512_SA", 512, 65_536),
            ("AAE_512_Transformer", "AAE_512", 512, 65_536),
            ("AAE_512+16", "AAE_512+16", 528, 67_584),
        ):
            entry = scores["models"][artifact_id]
            self.assertEqual(entry["model_id"], artifact_id)
            self.assertEqual(entry["display_name"], display_name)
            self.assertEqual(seal["models"][artifact_id]["display_name"], display_name)
            self.assertEqual(entry["latent"], {
                "tokens": tokens, "scalars": scalars, "decoder_tokens": 512, "decoder_scalars": 65_536,
            })
        for model, historical in historical_scores["models"].items():
            self.assertEqual(scores["models"][model]["sources"], historical["sources"])
        self.assertEqual(scores, self.cli("reload", *args, "--expected", score_path, script="compare.py"))
        self.assert_comparison_png(
            "figure", "--expected", score_path, self.root / "expanded-comparison.png",
        )

        seal["models"]["AAE_512+16"]["validation"]["sources"]["tenv"]["mse"] += 1
        corrupt_seal = self.root / "corrupt-seal.json"
        corrupt_seal.write_text(json.dumps(seal))
        rejected_scores = self.root / "rejected-scores.json"
        self.assertIn("selection seal verification failed", self.cli(
            "score", "--models", models, "--seal", corrupt_seal, "--output", rejected_scores,
            success=False, script="compare.py",
        ))
        self.assertFalse(rejected_scores.exists())

    def assert_invalid_expanded_selection(self, models, seal_path, seal_args):
        variant = "AAE_512+16"
        manifest_path = models / variant / "manifest.json"
        report_path = self.root / f"{variant}-validation.json"
        original_manifest = manifest_path.read_bytes()
        original_report = report_path.read_bytes()
        for field, replacement, error in (
            ("architecture", {}, "invalid AAE artifact"),
            ("recipe", {}, "invalid AAE artifact"),
            ("selected_epoch", 0, "checkpoint differs from earliest tenv minimum"),
        ):
            with self.subTest(field=field):
                manifest = json.loads(original_manifest)
                manifest["metadata"][field] = replacement
                manifest_path.write_text(json.dumps(manifest))
                report = json.loads(original_report)
                report["artifact_manifest_sha256"] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
                report_path.write_text(json.dumps(report))
                self.assertIn(error, self.cli("seal", *seal_args, success=False, script="compare.py"))
                self.assertFalse(seal_path.exists())
        manifest_path.write_bytes(original_manifest)
        for field, replacement, error in (
            ("reference_manifest_sha256", "wrong-mean", "does not match the artifact, reference and split"),
            ("provenance", {}, "does not match the artifact, reference and split"),
        ):
            with self.subTest(field=field):
                report = json.loads(original_report)
                report[field] = replacement
                report_path.write_text(json.dumps(report))
                self.assertIn(error, self.cli("seal", *seal_args, success=False, script="compare.py"))
                self.assertFalse(seal_path.exists())
        report = json.loads(original_report)
        report["sources"]["tenv"]["mse"] += 1
        report_path.write_text(json.dumps(report))
        self.assertIn("validation MSE differs", self.cli("seal", *seal_args, success=False, script="compare.py"))
        self.assertFalse(seal_path.exists())
        report_path.write_bytes(original_report)

    def test_expanded_figure_labels_each_attention_autoencoder(self):
        from matplotlib.figure import Figure

        from compare import plot_comparison

        metrics = {"profiles": 1, "pcc": 0.5, "rmse": 0.1, "relative_squared_error": 0.8}
        models = {}
        for artifact_id, display_name, tokens, scalars in (
            ("mean", "mean", None, None),
            ("PCA_32", "PCA_32", None, 32),
            ("PCA_128", "PCA_128", None, 128),
            ("PCA_512", "PCA_512", None, 512),
            ("MLP_32", "MLP_32", None, 32),
            ("MLP_128", "MLP_128", None, 128),
            ("MLP_512", "MLP_512", None, 512),
            ("AAE_512", "AAE_512_SA", 512, 65_536),
            ("AAE_512_Transformer", "AAE_512", 512, 65_536),
            ("AAE_512+16", "AAE_512+16", 528, 67_584),
            ("MatrixAE_64x8", "MatrixAE_64x8", 64, 512),
        ):
            latent = {"tokens": tokens, "scalars": scalars}
            if artifact_id.startswith("AAE_"):
                latent.update(decoder_tokens=512, decoder_scalars=65_536)
            models[artifact_id] = {
                "display_name": display_name, "latent": latent,
                "sources": {source: metrics for source in ("tahoe", "tenv")},
            }
        report = {"schema_version": 2, "role": "test", "models": models}
        with patch.object(Figure, "savefig", autospec=True) as savefig:
            plot_comparison(report, self.root / "comparison.png")
        figure = savefig.call_args.args[0]
        for subfigure in figure.subfigs:
            labels = [tick.get_text() for tick in subfigure.axes[0].get_yticklabels()]
            self.assertEqual(len(labels), 11)
            self.assertEqual([label for label in labels if label.startswith("AAE_")], [
                "AAE_512_SA  ·  512 tokens · 65,536 scalars",
                "AAE_512  ·  512 tokens · 65,536 scalars",
                "AAE_512+16  ·  528 tokens · 67,584 scalars",
            ])

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

    def test_matrix_ae_encodes_a_matrix_and_decodes_the_full_gene_axis(self):
        import torch
        from matrix_ae import MatrixAutoencoder

        model = MatrixAutoencoder(gene_count=3).eval()
        profiles = torch.tensor([[8.0, 3.0, 0.0], [0.0, 6.0, 4.0]])
        with torch.inference_mode():
            latent = model.encode(profiles)
            self.assertEqual(latent.shape, (2, 64, 8))
            reconstruction = model.decode(latent.clone())
            self.assertEqual(reconstruction.shape, profiles.shape)
            self.assertTrue(torch.isfinite(reconstruction).all())
            torch.testing.assert_close(model(profiles), reconstruction)
            torch.testing.assert_close(model.encode(profiles), latent, rtol=0, atol=0)

    def test_matrix_ae_cli_reloads_the_earliest_tenv_validation_minimum(self):
        split = self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        model_id = "MatrixAE_64x8"
        artifact = self.root / model_id
        fitted = self.cli("fit", "--model", model_id, "--output", artifact)
        metadata = fitted["metadata"]
        self.assertEqual(fitted["model_id"], model_id)
        self.assertEqual(metadata["latent_tokens"], 64)
        self.assertEqual(metadata["latent_scalars"], 512)
        self.assertEqual(metadata["architecture"]["latent_channels"], 8)
        self.assertEqual(metadata["fit_profiles"], {
            source: split["sources"][source]["roles"]["train"]["profiles"]
            for source in self.expressions
        })
        history = json.loads((artifact / "epochs.json").read_text())
        self.assertLessEqual(history[-1]["optimizer_steps"], 5)
        selected = min(history, key=lambda check: (check["validation_mse"]["tenv"], check["epoch"]))
        self.assertEqual(metadata["selected_epoch"], selected["epoch"])
        expected = self.root / f"{model_id}-validation.json"
        args = ["--artifact", artifact, "--mean-artifact", self.artifact, "--role", "validation"]
        scores = self.cli("score", *args, "--output", expected)
        self.assertEqual(scores["model_id"], model_id)
        self.assertEqual(scores["artifact_manifest_sha256"], fitted["artifact_manifest_sha256"])
        self.assertEqual(scores["provenance"], fitted["provenance"])
        for source, result in scores["sources"].items():
            self.assertEqual(result["profiles"], split["sources"][source]["roles"]["validation"]["profiles"])
            self.assertAlmostEqual(result["mse"], selected["validation_mse"][source])
        self.assertEqual(scores, self.cli("reload", *args, "--expected", expected))
        self.assert_comparison_png(
            "matrix-ae-loss-figure", "--history", artifact / "epochs.json", self.root / "matrix_ae_loss.png",
        )

    def test_matrix_mlp_encodes_a_matrix_and_decodes_the_full_gene_axis(self):
        import torch
        from matrix_mlp import MatrixMLPAutoencoder

        model = MatrixMLPAutoencoder(gene_count=3).eval()
        profiles = torch.tensor([[8.0, 3.0, 0.0], [0.0, 6.0, 4.0]])
        with torch.inference_mode():
            latent = model.encode(profiles)
            self.assertEqual(latent.shape, (2, 512, 128))
            reconstruction = model.decode(latent.clone())
            self.assertEqual(reconstruction.shape, profiles.shape)
            self.assertTrue(torch.isfinite(reconstruction).all())
            torch.testing.assert_close(model(profiles), reconstruction)
            torch.testing.assert_close(model(profiles[:1]), reconstruction[:1])

    def test_matrix_mlp_cli_reloads_the_earliest_tenv_validation_minimum(self):
        split = self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        model_id = "MLPm_512"
        artifact = self.root / model_id
        fitted = self.cli("fit", "--model", model_id, "--output", artifact)
        metadata = fitted["metadata"]
        self.assertEqual(fitted["model_id"], model_id)
        self.assertEqual(metadata["latent_tokens"], 512)
        self.assertEqual(metadata["latent_scalars"], 65_536)
        self.assertEqual(metadata["architecture"]["width"], 128)
        self.assertEqual(metadata["parameter_count"], 133_252)
        self.assertEqual(metadata["fit_profiles"], {
            source: split["sources"][source]["roles"]["train"]["profiles"]
            for source in self.expressions
        })
        history = json.loads((artifact / "epochs.json").read_text())
        self.assertEqual(history[-1]["optimizer_steps"], 5)
        selected = min(history, key=lambda check: (check["validation_mse"]["tenv"], check["epoch"]))
        self.assertEqual(metadata["selected_epoch"], selected["epoch"])
        expected = self.root / f"{model_id}-validation.json"
        args = ["--artifact", artifact, "--mean-artifact", self.artifact, "--role", "validation"]
        scores = self.cli("score", *args, "--output", expected)
        self.assertEqual(scores["artifact_manifest_sha256"], fitted["artifact_manifest_sha256"])
        self.assertEqual(scores["provenance"], fitted["provenance"])
        for source, result in scores["sources"].items():
            self.assertEqual(result["profiles"], split["sources"][source]["roles"]["validation"]["profiles"])
            self.assertAlmostEqual(result["mse"], selected["validation_mse"][source])
        self.assertEqual(scores, self.cli("reload", *args, "--expected", expected))

    def test_transformer_aae_encodes_reconstruction_tokens_and_decodes_the_full_gene_axis(self):
        import torch
        from aae import TransformerTranscriptomeAutoencoder

        model = TransformerTranscriptomeAutoencoder(gene_count=3).eval()
        profiles = torch.tensor([[8.0, 3.0, 0.0], [0.0, 6.0, 4.0]])
        with torch.inference_mode():
            latent = model.encode(profiles)
            self.assertEqual(latent.shape, (2, 512, 128))
            reconstruction = model.decode(latent)
            self.assertEqual(reconstruction.shape, profiles.shape)
            self.assertTrue(torch.isfinite(reconstruction).all())
            torch.testing.assert_close(model(profiles), reconstruction)

    def test_nuisance_queries_affect_reconstruction_through_the_encoder(self):
        import torch
        from aae import TransformerTranscriptomeAutoencoder

        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(20261109)
            model = TransformerTranscriptomeAutoencoder(gene_count=3, nuisance_tokens=16).eval()
        profiles = torch.tensor([[8.0, 3.0, 0.0], [0.0, 6.0, 4.0]])
        latent = model.encode(profiles)
        self.assertEqual(latent.shape, (2, 528, 128))
        reconstruction = model.decode(latent)
        self.assertEqual(reconstruction.shape, profiles.shape)
        torch.testing.assert_close(model(profiles), reconstruction)

        query_gradient, = torch.autograd.grad(reconstruction.square().sum(), model.latent_queries)
        self.assertTrue(torch.isfinite(query_gradient).all())
        self.assertTrue((query_gradient[:, 512:].norm(dim=-1) > 0).all())

        changed_latent = latent.detach().clone()
        changed_latent[:, 512:] = -changed_latent[:, 512:]
        torch.testing.assert_close(model.decode(changed_latent), reconstruction, rtol=0, atol=0)

    def test_transformer_aae_cli_reloads_the_earliest_tenv_validation_minimum(self):
        split = self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        recipes = []
        for model_id, tokens, scalars in (
            ("AAE_512_Transformer", 512, 65_536),
            ("AAE_512+16", 528, 67_584),
        ):
            artifact = self.root / model_id
            fitted = self.cli("fit", "--model", model_id, "--output", artifact)
            metadata = fitted["metadata"]
            self.assertEqual(fitted["model_id"], model_id)
            self.assertEqual(metadata["latent_tokens"], tokens)
            self.assertEqual(metadata["latent_scalars"], scalars)
            self.assertEqual(metadata["architecture"]["latent_tokens"], tokens)
            self.assertEqual(metadata["architecture"]["decoder_tokens"], 512)
            self.assertEqual(metadata["architecture"]["encoder_blocks"], 2)
            self.assertEqual(metadata["fit_profiles"], {
                source: split["sources"][source]["roles"]["train"]["profiles"]
                for source in self.expressions
            })
            recipes.append(metadata["recipe"])
            history = json.loads((artifact / "epochs.json").read_text())
            self.assertEqual(history[-1]["optimizer_steps"], 5)
            selected = min(history, key=lambda check: (check["validation_mse"]["tenv"], check["epoch"]))
            self.assertEqual(metadata["selected_epoch"], selected["epoch"])
            expected = self.root / f"{model_id}-validation.json"
            args = ["--artifact", artifact, "--mean-artifact", self.artifact, "--role", "validation"]
            scores = self.cli("score", *args, "--output", expected)
            self.assertEqual(scores["artifact_manifest_sha256"], fitted["artifact_manifest_sha256"])
            self.assertEqual(scores["provenance"], fitted["provenance"])
            for source, result in scores["sources"].items():
                expected_profiles = split["sources"][source]["roles"]["validation"]["profiles"]
                self.assertEqual(result["profiles"], expected_profiles)
                self.assertAlmostEqual(result["mse"], selected["validation_mse"][source])
            self.assertEqual(scores, self.cli("reload", *args, "--expected", expected))
        self.assertEqual(recipes[0]["objective"], "full-axis reconstruction MSE")
        self.assertEqual(recipes[0], recipes[1])

    def test_transformer_aae_rejects_incompatible_checkpoint_metadata(self):
        self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        model_id = "AAE_512_Transformer"
        artifact = self.root / model_id
        self.cli("fit", "--model", model_id, "--output", artifact)
        args = ["--artifact", artifact, "--mean-artifact", self.artifact]
        manifest_path = artifact / "manifest.json"
        manifest_bytes = manifest_path.read_bytes()
        for field, value, error in (
            ("architecture", {"gene_count": 3, "latent_tokens": 512, "width": 128, "heads": 8},
             "invalid AAE artifact"),
            ("recipe", {}, "invalid AAE artifact"),
            ("selected_epoch", 0, "differs from earliest tenv minimum"),
        ):
            with self.subTest(field=field):
                manifest = json.loads(manifest_bytes)
                manifest["metadata"][field] = value
                manifest_path.write_text(json.dumps(manifest))
                self.assertIn(error, self.cli("score", *args, success=False))

        manifest = json.loads(manifest_bytes)
        history_path = artifact / "epochs.json"
        history = json.loads(history_path.read_bytes())
        selected = history[manifest["metadata"]["selected_epoch"] - 1]
        selected["validation_mse"]["tahoe"] += 1
        manifest["metadata"]["selected_validation_mse"] = selected["validation_mse"]
        history_path.write_text(json.dumps(history))
        manifest["files"][history_path.name] = hashlib.sha256(history_path.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest))
        self.assertIn("validation MSE differs from the selected epoch",
                      self.cli("score", *args, success=False))

    def test_nuisance_aae_rejects_512_token_metadata_and_weights(self):
        self.cli("freeze-split")
        self.cli("fit", "--output", self.artifact)
        artifact = self.root / "AAE_512+16"
        self.cli("fit", "--model", "AAE_512+16", "--output", artifact)
        control = self.root / "AAE_512_Transformer"
        fitted_control = self.cli("fit", "--model", "AAE_512_Transformer", "--output", control)
        args = ["--artifact", artifact, "--mean-artifact", self.artifact]
        manifest_path = artifact / "manifest.json"
        manifest_bytes = manifest_path.read_bytes()
        for field in ("architecture", "latent_tokens", "latent_scalars"):
            with self.subTest(field=field):
                manifest = json.loads(manifest_bytes)
                manifest["metadata"][field] = fitted_control["metadata"][field]
                manifest_path.write_text(json.dumps(manifest))
                self.assertIn("invalid AAE artifact", self.cli("score", *args, success=False))

        manifest = json.loads(manifest_bytes)
        weights = (control / "weights.pt").read_bytes()
        (artifact / "weights.pt").write_bytes(weights)
        manifest["files"]["weights.pt"] = hashlib.sha256(weights).hexdigest()
        manifest_path.write_text(json.dumps(manifest))
        self.assertIn("size mismatch for latent_queries", self.cli("score", *args, success=False))

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
        self.assertEqual(history[-1]["optimizer_steps"], metadata["recipe"]["max_optimizer_steps"])
        self.assertEqual(metadata["stop_reason"], "optimizer step limit")
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
        self.assert_comparison_png(
            "aae-loss-figure", "--history", artifact / "epochs.json", self.root / "aae_loss.png",
        )

    def test_aae_stops_after_consecutive_non_improving_validation_checks(self):
        from aae import fit_aae, recipe_for
        from benchmark import Corpus, Split

        self.cli("freeze-split")
        corpus = Corpus(self.release, fixture=True)
        split = Split(corpus, self.split_path)
        artifact = self.root / "AAE_512"
        recipe = {**recipe_for(True, "cpu"), "max_optimizer_steps": 100}
        # A new best resets the count. Tied checkpoints keep the earlier epoch.
        tenv_losses = (4, 5, 5, 3, 3, 4, 5)
        checks = [{"tahoe": 1.0, "tenv": loss} for loss in tenv_losses]
        with patch("aae.recipe_for", return_value=recipe), \
                patch("neural.validation_mse", side_effect=checks) as validate:
            result = fit_aae(corpus, split, artifact, batch_size=None, device="cpu")
        self.assertEqual(validate.call_count, len(checks))
        self.assertEqual(result["metadata"]["stop_reason"], "validation patience")
        self.assertEqual(result["metadata"]["selected_epoch"], 4)

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
