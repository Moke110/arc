"""G1-G6 executable-gate tests for the v2 corpus (ticket #17).

Run from the repository root:

    python -m unittest discover -s stage2/02_dataset_corpus/execution/tests -v

The default suite is cheap: it checks the committed ``verify.json`` record and
the rendered README against the frozen v1 manifest, exercises the stage2 path
audit, and re-runs the disk-only structure/counts/crosswalk gates (G1-G3, a few
seconds; skipped when the local-only corpora are absent).  The two heavy gates
(G4 rebuild ~17 s / ~4 GiB, G5 double build ~60 s / ~8 GiB) are opt-in:

    ARC_RUN_SLOW_GATES=1 python -m unittest ... -v
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from cmdata import paths, verify  # noqa: E402

VERIFY_RECORD = paths.VERIFY_RECORD
RUN_SLOW = os.environ.get("ARC_RUN_SLOW_GATES") == "1"


def _missing_corpora() -> list[Path]:
    return [path for path in
            (verify._corpus_path(corpus_id)
             for corpus_id in verify.ALL_CORPORA) if not path.exists()]


def _corpora_skip() -> str:
    missing = _missing_corpora()
    return "" if not missing else f"corpus h5ads not present locally: {missing}"


def _frozen_skip() -> str:
    required = [verify.FROZEN_MANIFEST, verify.FROZEN_GENES,
                verify.FROZEN_VIEWS_PY, verify.FROZEN_PATHS_PY]
    missing = [path for path in required if not path.exists()]
    if missing:
        return f"frozen v1 evidence not present: {missing}"
    missing = [verify.FROZEN_VIEW_DIR / name for name in (
        "dtoxs_cm.h5ad", "screen_ctox.h5ad", "temp_env_atc-moa.h5ad",
        "temp_env_atc-moa-mesh.h5ad", "temp_env_all-smiles.h5ad")
        if not (verify.FROZEN_VIEW_DIR / name).exists()]
    return "" if not missing else f"frozen v1 view h5ads not present: {missing}"


def _record_skip() -> str:
    return "" if VERIFY_RECORD.exists() else (
        f"verify record not written yet: {VERIFY_RECORD}")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _frozen_manifest() -> dict[str, dict[str, str]]:
    with open(verify.FROZEN_MANIFEST, newline="") as handle:
        return {row["view"]: row for row in csv.DictReader(handle)}


@unittest.skipUnless(not _record_skip(), _record_skip())
class VerifyRecordTests(unittest.TestCase):
    """The committed verify.json + README gate block are internally consistent."""

    @classmethod
    def setUpClass(cls):
        cls.record = json.loads(VERIFY_RECORD.read_text())
        cls.manifest = _frozen_manifest()

    def test_all_six_gates_pass(self):
        self.assertEqual(self.record["ticket"], "#17")
        self.assertEqual(self.record["result"], "PASS")
        self.assertEqual(self.record["total"], 6)
        self.assertEqual(self.record["passed"], 6)
        self.assertEqual(set(self.record["gates"]),
                         set(verify.GATE_ORDER))
        for gate in verify.GATE_ORDER:
            self.assertTrue(self.record["gates"][gate]["passed"], gate)

    def test_g4_five_hashes_equal_the_frozen_manifest(self):
        g4 = self.record["gates"]["G4"]["data"]
        self.assertEqual(g4["matched"], 5)
        self.assertEqual(g4["expected"], 5)
        self.assertEqual(len(g4["views"]), 5)
        seen = set()
        for view in g4["views"]:
            self.assertTrue(view["match"], view["view"])
            self.assertRegex(view["rebuilt_content_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(view["rebuilt_content_sha256"],
                             view["frozen_content_sha256"])
            self.assertEqual(view["rebuilt_content_sha256"],
                             self.manifest[view["view"]]["content_sha256"])
            seen.add(view["view"])
        self.assertEqual(seen, set(self.manifest))

    def test_g5_records_a_double_build(self):
        data = self.record["gates"]["G5"]["data"]
        self.assertEqual(data["corpus"], "dtoxs_ipsc")
        self.assertEqual(data["runs"], 2)
        self.assertEqual(data["file_sha256_run1"],
                         data["file_sha256_run2"])
        self.assertEqual(data["content_sha256"],
                         data["record_content_sha256"])
        self.assertTrue(all(data["checks"].values()))

    def test_g6_records_zero_stale(self):
        data = self.record["gates"]["G6"]["data"]
        self.assertEqual(data["exit_code"], 0)
        self.assertEqual(data["stale"], 0)
        self.assertGreater(data["files_scanned"], 0)

    def test_readme_gate_section_is_rendered_from_the_record(self):
        text = paths.CORPUS_V2_README.read_text()
        block = "\n".join(verify.render_gate_block(self.record))
        self.assertIn(block, text,
                      "README gate section is stale; re-run --verify")
        for gate in verify.GATE_ORDER:
            self.assertIn(verify.GATE_LABELS[gate], text)
        for view in self.record["gates"]["G4"]["data"]["views"]:
            self.assertIn(view["rebuilt_content_sha256"], text)

    def test_coverage_record_matches_the_readme(self):
        coverage = json.loads(paths.COVERAGE_RECORD.read_text())
        self.assertEqual(coverage["readme"]["sha256"],
                         _sha256_file(paths.CORPUS_V2_README))
        self.assertEqual(coverage["gate_results"]["verified_by"],
                         verify._rel(VERIFY_RECORD))
        self.assertEqual(coverage["verification"]["result"], "PASS")


class PathAuditTests(unittest.TestCase):
    """G6's script reports 0 stale literals on the tracked tree."""

    def test_check_stage2_paths_is_clean(self):
        script = paths.REPO_ROOT / "script" / "check_stage2_paths.sh"
        process = subprocess.run(["bash", str(script)], cwd=str(paths.REPO_ROOT),
                                 capture_output=True, text=True, timeout=300)
        self.assertEqual(process.returncode, 0, process.stdout
                         + process.stderr)
        self.assertIn("stale: 0 occurrences", process.stdout)
        self.assertIn("RESULT: PASS", process.stdout)


@unittest.skipUnless(not _corpora_skip(), _corpora_skip())
class FastGateTests(unittest.TestCase):
    """G1-G3 re-derived from the corpus files on disk."""

    def test_g1_structure(self):
        check = verify.check_g1_structure()
        self.assertTrue(check.passed, check.data["problems"])
        self.assertEqual(set(check.data["corpora"]),
                         set(verify.ALL_CORPORA))

    def test_g2_counts(self):
        check = verify.check_g2_counts()
        self.assertTrue(check.passed, check.data["problems"])
        self.assertEqual(set(check.data["corpora"]),
                         set(verify.ALL_CORPORA))

    def test_g3_crosswalk(self):
        check = verify.check_g3_crosswalk()
        self.assertTrue(check.passed, check.data["problems"])
        self.assertEqual(check.data["totals"]["genes_symbols"], 18164)
        for corpus_id in (*paths.CORPUS_IDS, "tahoe"):
            self.assertIn(corpus_id, check.data["totals"])


@unittest.skipUnless(not _corpora_skip() and not _frozen_skip(),
                     _corpora_skip() or _frozen_skip())
@unittest.skipUnless(RUN_SLOW, "set ARC_RUN_SLOW_GATES=1 for G4/G5")
class SlowGateTests(unittest.TestCase):
    """The heavy gates: G4's v1 rebuild and G5's double corpus build."""

    def test_g4_reconstructs_the_five_frozen_views(self):
        check = verify.check_g4_v1_reconstruction()
        self.assertTrue(check.passed, check.data["problems"])
        self.assertEqual(check.data["matched"], 5)
        manifest = _frozen_manifest()
        for view in check.data["views"]:
            self.assertEqual(view["rebuilt_content_sha256"],
                             manifest[view["view"]]["content_sha256"])

    def test_g5_rebuild_is_byte_reproducible(self):
        check = verify.check_g5_reproducibility()
        self.assertTrue(check.passed, check.data["problems"])
        self.assertTrue(all(check.data["checks"].values()))


if __name__ == "__main__":
    unittest.main()
