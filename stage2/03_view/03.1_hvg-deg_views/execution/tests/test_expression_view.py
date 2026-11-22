"""Worked examples through published view builder and reader."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

STEP = Path(__file__).resolve().parents[2]
STAGE = STEP.parents[1]
sys.path[:0] = [str(item) for item in (
    STEP / "execution", STEP / "output/api",
    STAGE / "01_corpus/output/api",
    STAGE / "02_axis_view/02.1_hvg-deg/output/api",
)]

from native_corpus import open_corpus, seal_corpus
from prediction_axis import COMMON_SHA256, PredictionAxis
from build_views import publish_views, replay_views
from expression_view import open_view


class ExpressionViewTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        corpus = self.root / "counts"
        corpus.mkdir()
        roles = ["treated", "vehicle", "vehicle", "treated", "treated", "vehicle", "treated"]
        obs = pd.DataFrame({
            "role": roles, "matched_control_group": ["A", "A", "A", "A", "B", "B", ""],
            "cell_line": "CM", "drug_parent_inchikey": ["P" if r == "treated" else "" for r in roles],
            "dose_um": [1.0 if r == "treated" else 0.0 for r in roles],
            "view_exclusion_reason": "",
        }, index=pd.Index(["t1", "v1", "v2", "t2", "t3", "v3", "unmatched"], name="profile_id"))
        vectors = np.zeros((7, 18164), dtype=np.int32)
        vectors[:, 1] = 100000
        vectors[:, [5, 2]] = [[3, 0], [1, 0], [0, 9], [0, 2], [1, 4], [10, 0], [1, 1]]
        for source in ("tenv", "tahoe_bulk"):
            AnnData(X=vectors, obs=obs, var=pd.DataFrame(index=[str(i) for i in range(18164)])).write_h5ad(corpus / f"{source}.h5ad")
        seal_corpus(corpus, dict.fromkeys(("tenv", "tahoe_bulk"), dict(enumerate(range(18164)))),
                    18164, COMMON_SHA256, {}, {})
        self.corpus = self.enterContext(open_corpus(corpus))
        self.axis = PredictionAxis((5, 2), ("five", "two"), "a" * 64)
        self.output = self.root / "views"

    def test_raw_first_control_and_individual_vehicle_rows(self):
        publish_views(self.output, self.corpus, self.axis)
        with open_view(self.output, source="tenv") as view:
            self.assertEqual(list(view.obs.index), ["t1", "v1", "v2", "t2", "t3", "v3"])
            self.assertEqual(list(view.var.gene_index), [5, 2])
            np.testing.assert_allclose(view.read(0, 1, "control")[0], np.log1p([1000, 9000]))
            np.testing.assert_allclose(view.read(1, 3), np.log1p([[10000, 0], [0, 10000]]))
            self.assertEqual(view.exclusions.set_index("profile_id").loc["unmatched", "reason"],
                             "no_matched_vehicle")
        with open_view(self.output, source="tenv", family="merge") as view:
            np.testing.assert_allclose(view.read(0, 1, "control")[0], np.log1p([5500, 4500]))
            np.testing.assert_allclose(view.read(0, 1)[0], np.log1p([4000, 6000]))
            np.testing.assert_allclose(view.read(0, 1, "delta"), view.read(0, 1) - view.read(0, 1, "control"))

    def test_families_and_pooled_label_release_are_replayable(self):
        publish_views(self.output, self.corpus, self.axis)
        with open_view(self.output, source="tenv", family="full-pooled-label") as view:
            np.testing.assert_allclose(view.read(0, 1)[0], np.log1p([10000, 0]))
            np.testing.assert_allclose(view.read(0, 1, "control")[0], np.log1p([1000, 9000]))
            delta = np.log1p([4000, 6000]) - np.log1p([5500, 4500])
            np.testing.assert_allclose(view.read(0, 1, "delta")[0], delta, atol=1e-6)
            np.testing.assert_allclose(view.read(3, 5, "delta"), [delta, delta], atol=1e-6)
        with open_view(self.output, source="tenv", family="solo-deepest") as view:
            self.assertEqual(list(view.obs.index), ["t3"])
            np.testing.assert_allclose(view.read(0, 1, "control")[0], np.log1p([10000, 0]))
        for source in ("tenv", "tahoe_bulk"):
            with open_view(self.output, source=source, family="solo") as view:
                self.assertEqual(view.shape, (1, 2))
                self.assertIn(view.obs.index[0], {"t1", "t2", "t3"})
        with self.assertRaisesRegex(ValueError, "unknown expression view arm"):
            open_view(self.output, source="tahoe_bulk", family="solo-deepest")
        self.assertEqual(replay_views(self.output, self.corpus, self.axis)["result"], "PASS")

    def test_reader_rejects_stale_digest_and_modified_payload(self):
        publish_views(self.output, self.corpus, self.axis)
        with self.assertRaisesRegex(ValueError, "generation"):
            open_view(self.output, source="tenv", corpus_sha256="0" * 64)
        with open_view(self.output, source="tenv") as view:
            with self.assertRaisesRegex(ValueError, "expression window"):
                view.read(0, 257)
        path = self.output / "full/tenv.h5ad"
        with path.open("ab") as output:
            output.write(b"unexpected")
        with self.assertRaisesRegex(ValueError, "payload digest"):
            open_view(self.output, source="tenv")
