import unittest

import numpy as np

from tsg_drug_emb.representations import checked_vectors


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
