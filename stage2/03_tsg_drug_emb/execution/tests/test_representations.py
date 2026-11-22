import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

from tsg_drug_emb.representations import (
    ARMS, PARENT_COLUMNS, SCHEMA, WIDTHS, RepresentationInputs, checked_vectors, fold_coverage, sha256,
)


class RepresentationContract(unittest.TestCase):
    def test_native_values_and_width_are_preserved(self):
        native = np.full((2, 2304), -12.75, dtype=np.float32)
        self.assertIs(checked_vectors("KPGT", native, 2), native)

    def test_incomplete_nonfinite_or_nonfloat32_inputs_are_rejected(self):
        for bad in (
            np.ones((1, 2304), dtype=np.float32),
            np.ones((2, 2304), dtype=np.float64),
            np.full((2, 2304), np.nan, dtype=np.float32),
        ):
            with self.subTest(shape=bad.shape, dtype=bad.dtype):
                with self.assertRaisesRegex(ValueError, "KPGT"):
                    checked_vectors("KPGT", bad, 2)

    def test_fingerprints_reject_nonbinary_values(self):
        invalid = np.full((2, 1024), 0.5, dtype=np.float32)
        with self.assertRaisesRegex(ValueError, "binary"):
            checked_vectors("E3FP", invalid, 2)


class TahoeRepresentationContract(unittest.TestCase):
    def test_sealed_tahoe_vectors_are_read_in_frozen_parent_order(self):
        keys = ["A", "B", "C"]
        parents = pd.DataFrame(zip(keys, ("C", "N", "O")), columns=PARENT_COLUMNS)
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "tahoe_bulk"
            source.mkdir()
            (source / "manifest.json").write_text("{}")
            inputs = SimpleNamespace(
                source="tahoe_bulk", directory=source, obs=parents, parents=parents,
                parent_lookup=dict(zip(keys, range(3))),
                rows=lambda scenario, fold, partition: np.array(
                    [(fold + ("test", "train", "validation").index(partition)) % 3]),
                morgan=np.zeros((3, WIDTHS["ECFP4C"]), dtype=np.float32),
            )
            freeze = root / "features"
            freeze.mkdir()
            parents.to_csv(freeze / "parents.csv", index=False)
            parents[[PARENT_COLUMNS[0]]].assign(conformers=1).to_csv(
                freeze / "e3fp_conformers.csv", index=False)
            for arm in ARMS:
                values = np.zeros((3, WIDTHS[arm]), dtype=np.float32)
                if arm == "KPGT":
                    values[:, 0] = [-9, 4, 7]
                np.save(freeze / f"{arm}.npy", values)
            manifest = dict(schema=SCHEMA, source=inputs.source, widths=WIDTHS,
                            coverage=fold_coverage(inputs), inputs_sha256=sha256(source / "manifest.json"),
                            files={path.name: sha256(path) for path in freeze.iterdir()})
            (freeze / "manifest.json").write_text(json.dumps(manifest))

            freeze_inputs = RepresentationInputs(inputs, freeze)
            np.testing.assert_array_equal(freeze_inputs.vectors("KPGT", np.array([2, 0]))[:, 0],
                                          np.array([7, -9], dtype=np.float32))
