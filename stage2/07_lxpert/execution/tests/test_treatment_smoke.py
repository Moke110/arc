"""Exercise the treatment integration through the Step 07 CLI."""
import hashlib
import json
from math import log10
import shutil
import unittest

import numpy as np

import test_formal_fit


class TreatmentSmokeTest(unittest.TestCase):
    setUp = test_formal_fit.FormalFitTest.setUp
    cli = test_formal_fit.FormalFitTest.cli

    def test_treatment_smoke_uses_measured_inputs_and_frozen_checkpoint(self):
        self.cli("freeze", "--batch-size", "2", "--max-epochs", "1", "--max-seconds", "60")
        self.cli("develop")
        output = self.release.parent / "autoencoder"
        refit = self.cli("refit", "--output", str(output))
        tenv = self.release.parent / "matched"
        tenv.mkdir()
        shutil.copyfile(self.release / "genes.csv", tenv / "genes.csv")
        (tenv / "wells.csv").write_text(
            "case_id,drug_key,dose_um,control_index,group_id\n"
            "well-a,drug-a,0.5,1,group-a\n"
            "well-b,drug-b,0.05,0,group-b\n"
        )
        (tenv / "controls.csv").write_text("control_index,plate,n_wells\n0,p0,2\n1,p1,3\n")
        (tenv / "drugs.csv").write_text("drug_key,morgan_index\ndrug-b,0\ndrug-a,1\n")
        baseline = np.array([3, 1, 2], dtype=np.float32)
        np.save(tenv / "controls.npy", np.stack([np.zeros_like(baseline), baseline]))
        morgan = np.zeros((2, 2048), dtype=np.uint8)
        morgan[1, [0, 123, 2047]] = 1
        np.save(tenv / "morgan.npy", morgan)
        files = {file.name: hashlib.sha256(file.read_bytes()).hexdigest()
                 for file in tenv.iterdir()}
        (tenv / "manifest.json").write_text(json.dumps({
            "schema": "tenv-input-v1", "gene_count": 3, "drug_count": 2, "files": files,
            "morgan": {"radius": 2, "bits": 2048, "chirality": True, "binary": True,
                       "invariants": "atom"},
        }))

        args = ("--tenv-input", str(tenv), "--checkpoint", str(output / "checkpoint.pt"),
                "--case-id", "well-a")
        smoke = self.cli("smoke", *args)
        self.assertEqual(smoke["checkpoint_sha256"], refit["checkpoint_sha256"])
        self.assertEqual(smoke["fit_stage"], "full_data_refit")
        self.assertEqual(smoke["case"]["plate"], "p1")
        self.assertEqual(smoke["case"]["dmso_wells"], 3)
        self.assertEqual(smoke["case"]["dose_um"], 0.5)
        self.assertEqual(smoke["case"]["log10_relative_dose"], log10(10))
        self.assertEqual(smoke["case"]["morgan_active_bits"], 3)
        self.assertEqual(smoke["case"]["baseline_sha256"],
                         hashlib.sha256(baseline.tobytes()).hexdigest())
        self.assertEqual(smoke["latent_shape"], [1, 512, 128])
        self.assertEqual(smoke["condition_shape"], [1, 2, 128])
        self.assertEqual(smoke["block_order"], ["CA", "SA", "SA", "CA"])
        self.assertEqual(smoke["post_treatment_shape"], [1, 3])
        self.assertEqual(smoke["delta_shape"], [1, 3])
        self.assertEqual(smoke["delta_max_abs_error"], 0.0)
        self.assertTrue(smoke["finite"])
        self.assertTrue(smoke["autoencoder_unchanged"])
        self.assertTrue(smoke["autoencoder_gradients_absent"])
        for name in ("morgan_adapter", "dose_adapter"):
            self.assertGreater(smoke["gradient_l2"][name], 0.0)

        (tenv / "controls.npy").write_bytes(b"tampered")
        self.cli("smoke", *args, error="tenv input digest mismatch: controls.npy")
