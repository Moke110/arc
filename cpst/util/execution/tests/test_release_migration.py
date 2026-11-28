"""Exercise CPST's migrated receipt and historical source interfaces."""

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import unittest

from cpst_util.release_names import (checkpoint_registration_sha256, current_path, historical_path,
                                     original_axis_manifest_sha256)

ROOT = Path(__file__).resolve().parents[4]


class ReleaseNamesTest(unittest.TestCase):
    def test_round_trip_a_published_source(self):
        original = "stage2_CarPS-VAE/07_cpsv_alpha/output/api/carps_vae_alpha.py"
        current = "/".join(("cpst", "07_cpst_alpha", "output/api/cpst_alpha.py"))
        self.assertEqual(current_path(original), current)
        self.assertEqual(historical_path(current), original)
        self.assertFalse((ROOT / current).exists())
        subprocess.check_call(
            ["git", "cat-file", "-e", f"a7613b74fc23f9f48d8ea0d4479542c008d701db:{current}"], cwd=ROOT,
        )

    def test_historical_hash_uses_the_original_git_blob(self):
        from util.training.provenance import git_code_hashes

        revision = "825f48c96a0968c11b0bffa9a8b06686b0d16e03"
        current = "/".join(("cpst", "07_cpst_alpha", "output/api/cpst_alpha.py"))
        original = "stage2_CarPS-VAE/07_cpsv_alpha/output/api/carps_vae_alpha.py"
        frozen = subprocess.check_output(["git", "show", f"{revision}:{original}"], cwd=ROOT)
        self.assertEqual(git_code_hashes(ROOT, revision, (current,)), {current: sha256(frozen).hexdigest()})

    def test_historical_figure_receipt_reissues_the_renderers_and_preserves_plotted_bytes(self):
        old = "stage2_CarPS-VAE/09_cpsv_axis/log/figure-render.json"
        revision = "9e7c85e3"
        original = json.loads(subprocess.check_output(["git", "show", f"{revision}:{old}"], cwd=ROOT))
        path = current_path(old)
        current = json.loads((ROOT / path).read_text())
        evidence = current["evidence"]
        predecessor = evidence["predecessor_revision"]
        frozen = subprocess.check_output(["git", "show", f"{predecessor}:{path}"], cwd=ROOT)
        self.assertEqual(sha256(frozen).hexdigest(), evidence["predecessor_figure_receipt_sha256"])
        renamed = json.loads(frozen)
        self.assertEqual(renamed["schema"], "cpst-figure-render-v1")
        self.assertEqual(renamed["step"], "09_cpst_axis")
        for field in ("csv_files_sha256", "png_files_sha256", "environment"):
            self.assertEqual(renamed[field], original[field])
        renderers = renamed["renderer_files_sha256"]
        self.assertEqual(set(renderers), {current_path(path) for path in original["renderer_files_sha256"]})
        for path, digest in renderers.items():
            frozen_renderer = subprocess.check_output(["git", "show", f"{predecessor}:cpst/{path}"], cwd=ROOT)
            self.assertEqual(sha256(frozen_renderer).hexdigest(), digest)

    def test_migrated_checkpoint_keeps_its_embedded_registration_hash(self):
        migration = json.loads((ROOT / "cpst/util/log/385-release-migration.json").read_bytes())
        suffix = "/Zmol-FiLM/fold-0/seed-17/registration.json"
        record, = [row for row in migration["metadata"] if row["path"].endswith(suffix)]
        self.assertEqual(
            checkpoint_registration_sha256(record["sha256"]),
            "ec67859b8dc8c4538f40913ca4fbef50c4ab9e0f2c6e1da74649bd1891e173a3",
        )
        self.assertEqual(checkpoint_registration_sha256("new-registration"), "new-registration")

    def test_only_the_published_denoise_manifest_resolves_to_its_original_hash(self):
        relative = "02_gene_axis/02.2_denoise-deg/output/release/manifest.json"
        current = (ROOT / "cpst" / relative).read_bytes()
        original = subprocess.check_output(
            ["git", "show", f"9e7c85e3:stage2_CarPS-VAE/{relative}"], cwd=ROOT,
        )
        observed = sha256(current).hexdigest()
        self.assertEqual(original_axis_manifest_sha256(observed), sha256(original).hexdigest())
        tampered = sha256(current + b"\n").hexdigest()
        self.assertEqual(original_axis_manifest_sha256(tampered), tampered)

    def test_path_matching_does_not_rewrite_unrelated_names(self):
        unaffected = "fixtures/not-stage2_CarPS-VAE/archive_07_cpsv_alpha.py"
        self.assertEqual(current_path(unaffected), unaffected)
        self.assertEqual(historical_path("fixtures/not-cpst/archive_07_cpst_alpha.py"),
                         "fixtures/not-cpst/archive_07_cpst_alpha.py")


if __name__ == "__main__":
    unittest.main()
