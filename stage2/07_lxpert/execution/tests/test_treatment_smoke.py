"""Exercise the treatment integration through the Step 07 CLI."""
import hashlib
import json
from math import log10
import shutil
import subprocess
import sys
import unittest

import numpy as np

import test_formal_fit
import test_inputs


class TreatmentSmokeTest(unittest.TestCase):
    setUp = test_formal_fit.FormalFitTest.setUp
    cli = test_formal_fit.FormalFitTest.cli

    def test_smoke_default_checkpoint_requires_the_canonical_digest(self):
        repo = self.release.parent / "checkout"
        step = repo / "stage2/07_lxpert"
        (repo / ".git").mkdir(parents=True)
        scripts = step / "execution/script"
        shutil.copytree(test_inputs.CLI.parent, scripts)
        records = repo / "stage2/02_dataset_corpus/execution/build_record"
        records.mkdir(parents=True)
        (records / "tahoe_corpus_move.json").write_text('{"file_sha256": "unused"}')
        (records / "corpus_build.json").write_text('{"corpora": {"tenv": {"file_sha256": "unused"}}}')
        record = step / "execution/formal_fit_record/run.json"
        record.parent.mkdir(parents=True)
        record.write_text(json.dumps({"checkpoint_sha256": hashlib.sha256(b"published").hexdigest()}))
        checkpoint = step / "output/autoencoder/checkpoint.pt"
        checkpoint.parent.mkdir(parents=True)
        checkpoint.write_bytes(b"changed")

        for args in ((), ("--checkpoint", str(checkpoint.relative_to(repo)))):
            with self.subTest(args=args):
                result = subprocess.run(
                    [sys.executable, str(scripts / "lxpert.py"), "smoke", "--fixture",
                     "--device", "cpu", "--release", str(self.release), *args],
                    cwd=repo, capture_output=True, text=True, check=False,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("checkpoint differs from the verified Step 07 autoencoder release",
                              result.stderr)

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

        controls = tenv / "controls.npy"
        original_controls = controls.read_bytes()
        controls.write_bytes(b"tampered")
        self.cli("smoke", *args, error="tenv input digest mismatch: controls.npy")
        controls.write_bytes(original_controls)

        refit_record = output / "refit.json"
        wrong_mode = json.loads(refit_record.read_text())
        wrong_mode["earlier_tenv_test_expression_in_refit"] = True
        refit_record.write_text(json.dumps(wrong_mode))
        (output / "refit.sha256").write_text(hashlib.sha256(refit_record.read_bytes()).hexdigest())
        self.cli("smoke", *args, error="formal checkpoint provenance mismatch")

    def test_smoke_rejects_a_short_fit_checkpoint(self):
        output = self.release.parent / "short_fit"
        self.cli("fit", "--steps", "1", "--batch-size", "2", "--output", str(output))
        self.cli("smoke", "--checkpoint", str(output / "checkpoint.pt"),
                 error="treatment smoke requires a full-data refit checkpoint")
