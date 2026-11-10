"""Fixture validation and replay for the frozen expression interface and matrix PCA."""

import unittest
from unittest.mock import patch

import test_benchmark


class SharedExpressionTest(unittest.TestCase):
    def setUp(self):
        self.fixture = test_benchmark.BenchmarkTest()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.setUp()
        self.fixture.cli("freeze-split")
        self.fixture.cli("fit", "--output", self.fixture.artifact)

    def test_interface_validation_and_replay(self):
        fixture = self.fixture
        interface = fixture.root / "interface"
        record = fixture.cli("fit", "--model", "shared_expression", "--output", interface)
        self.assertEqual(record["metadata"]["fit_role"], "train")
        report = fixture.root / "interface-validation.json"
        options = ("--artifact", interface, "--mean-artifact", fixture.artifact)
        validation = fixture.cli("score", *options, "--output", report)
        self.assertEqual(validation, fixture.cli("reload", *options, "--expected", report))
        self.assertEqual(set(validation["sources"]), {"tahoe", "tenv"})
        self.assertIn("validation", fixture.cli("score", *options, "--role", "test",
                                               "--allow-test", success=False))

    def test_interface_calibration_reads_only_joint_training_membership(self):
        from benchmark import Corpus, Split
        from shared_expression import fit_interface

        fixture = self.fixture
        corpus = Corpus(fixture.release, fixture=True)
        split = Split(corpus, fixture.split_path)
        read = corpus.read
        visited = set()

        def train_only(source, rows):
            self.assertLessEqual(set(rows), set(split.rows[source]["train"]))
            visited.add(source)
            return read(source, rows)

        with patch.object(corpus, "read", side_effect=train_only):
            fit_interface(corpus, split, fixture.root / "interface", "cpu")
        self.assertEqual(visited, {"tahoe", "tenv"})

    def test_matrix_pca_reuses_one_basis_for_all_channels(self):
        from pathlib import Path

        import numpy as np
        import torch

        from benchmark import Corpus, Split
        from shared_expression import load_interface
        from shared_pca import MatrixPCA

        fixture = self.fixture
        interface = fixture.root / "interface"
        fixture.cli("fit", "--model", "shared_expression", "--output", interface)
        artifact = fixture.root / "pca"
        record = fixture.cli("fit", "--model", "PCA_1024", "--interface", interface,
                             "--output", artifact)
        self.assertEqual(record["metadata"]["latent_shape"], [3, 16])

        corpus = Corpus(fixture.release, fixture=True)
        split = Split(corpus, fixture.split_path)
        _, model = load_interface(interface, corpus, split, "cpu")
        train = np.concatenate([corpus.read(source, split.rows[source]["train"])
                                for source in ("tahoe", "tenv")])
        representation = model.generate(torch.from_numpy(train)).numpy()
        center = np.load(artifact / "center.npy")
        np.testing.assert_allclose(center, representation.mean(axis=0), atol=2e-7)
        components = np.load(artifact / "components.npy")
        np.testing.assert_allclose(components @ components.T, np.eye(3), atol=1e-6)
        pooled = (representation - center).transpose(0, 2, 1).reshape(-1, 3)
        covariance = pooled.T @ pooled
        covariance /= np.trace(covariance)
        eigenvalues = np.linalg.eigvalsh(covariance)[::-1]
        np.testing.assert_allclose(components @ covariance @ components.T,
                                   np.diag(eigenvalues), atol=1e-6)

        projection = MatrixPCA(torch.from_numpy(center), torch.from_numpy(components))
        value = model.generate(torch.tensor([[1., 2., 3.]]))
        encoded = projection.encode(value)
        np.testing.assert_allclose(encoded.numpy(), components @ (value.numpy() - center),
                                   atol=1e-6)
        np.testing.assert_allclose(projection.decode(encoded).numpy(), value.numpy(), atol=1e-6)
        options = ("--artifact", artifact, "--interface", interface, "--mean-artifact", fixture.artifact)
        report = fixture.root / "pca-validation.json"
        validation = fixture.cli("score", *options, "--output", report)
        self.assertEqual(validation, fixture.cli("reload", *options, "--expected", report))
        other = fixture.root / "other-interface"
        fixture.cli("fit", "--model", "shared_expression", "--output", other)
        manifest = Path(other) / "manifest.json"
        manifest.write_text(manifest.read_text() + "\n")
        self.assertIn("interface", fixture.cli("score", *options, "--interface", other, success=False))

        axis = fixture.release / "genes.csv"
        header, *genes = axis.read_text().splitlines()
        axis.write_text("\n".join([header, *reversed(genes)]) + "\n")
        fixture.write_manifest()
        fixture.split_path = fixture.root / "other-split.json"
        fixture.cli("freeze-split")
        self.assertIn("axis", fixture.cli("score", *options, success=False))
