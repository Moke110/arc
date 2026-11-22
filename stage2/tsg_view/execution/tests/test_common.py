"""Shared-implementation seams: canonical JSON, repo discovery, schedule and population.

These are the single implementations the step modules import. The tests pin the two
properties the rest of the step relies on: the canonical JSON form is byte-identical to
the published manifests, and the released-population constants reconcile to the counts
the three consumers (``views.merge``, ``views.full_pooled_label``, ``truth``) assert.
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import common  # noqa: E402
import matrix  # noqa: E402
import render  # noqa: E402
import truth  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
VIEWS = REPO / "stage2/tsg_view/output/views"


class CanonicalJsonTest(unittest.TestCase):
    def test_manifest_round_trip_is_byte_identical(self):
        """The one serializer must reproduce every published view manifest exactly.

        A de-duplicated helper that changed indentation, key order, the trailing newline
        or the number formatting would silently change a sealed product digest, so this
        runs the published manifests through the helper and compares raw bytes.
        """
        self.assertEqual(common.sha256.__module__, "common")
        manifests = sorted(VIEWS.glob("*/*.json"))
        if not manifests:
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "value.json"
                common.write_json(path, {"b": 1, "a": [2, 3]})
                self.assertEqual(path.read_text(), '{\n  "a": [\n    2,\n    3\n  ],\n'
                                                   '  "b": 1\n}\n')
                self.assertEqual(common.read_json(path), {"a": [2, 3], "b": 1})
                return
        with tempfile.TemporaryDirectory() as tmp:
            for path in manifests:
                written = Path(tmp) / path.name
                common.write_json(written, common.read_json(path))
                self.assertEqual(written.read_bytes(), path.read_bytes(), str(path))

    def test_non_finite_values_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                common.write_json(Path(tmp) / "bad.json", {"pcc": float("nan")})


class RepositoryRootTest(unittest.TestCase):
    def test_discovery_agrees_with_the_importing_module(self):
        self.assertEqual(common.repo_root(), REPO)
        self.assertEqual(Path(common.repo_root(REPO)), REPO)
        with self.assertRaisesRegex(ValueError, "repository root"):
            common.repo_root(HERE)

    def test_relative_falls_back_to_the_bare_name(self):
        self.assertEqual(common.relative(REPO / "some/dir/x.csv", REPO), "some/dir/x.csv")
        self.assertEqual(common.relative("/elsewhere/y.csv", REPO), "y.csv")


class ScheduleTruthTest(unittest.TestCase):
    def test_one_view_order_and_seed_source(self):
        """matrix, render and report must read the same schedule, not copies of it."""
        self.assertEqual(common.SOURCE_ORDER, matrix.SOURCE_ORDER)
        self.assertEqual(common.VIEW_ORDER, matrix.VIEW_ORDER)
        self.assertEqual(common.VIEW_ORDER, render.VIEW_ORDER)
        for source in common.SOURCE_ORDER:
            for view in common.VIEW_ORDER[source]:
                self.assertEqual(common.seeds_for(view), matrix.seeds_for(view))
                self.assertEqual(common.seeds_for(view), render.SEEDS[source][view])
        self.assertEqual(common.seeds_for("full"), (42, 43))
        self.assertEqual(common.seeds_for("solo-deepest"), (42,))
        self.assertEqual(common.SEEDS, matrix.SEEDS)


class ReleasedPopulationTest(unittest.TestCase):
    """The published solo/solo-deepest member counts are the released histogram."""

    def test_consumers_read_the_same_constants(self):
        import views.merge as merge
        import views.full_pooled_label as pooled

        for source in common.SOURCE_ORDER:
            released = common.RELEASED[source]
            self.assertEqual(common.released_conditions()[source], released["conditions"])
            self.assertEqual(common.released_histogram(source), released["member_histogram"])
            self.assertEqual(sum(released["member_histogram"].values()),
                             released["conditions"])
            self.assertEqual(common.replicated_conditions()[source],
                             released["conditions"] - released["member_histogram"]["1"])
            self.assertEqual(merge.EXPECTED_CONDITIONS[source], released["conditions"])
            self.assertEqual({key: value for key,
                              value in common.released_histogram(source).items() if value},
                             merge.EXPECTED_MEMBER_HISTOGRAM[source])
            self.assertEqual(pooled.EXPECTED_CONDITIONS[source], released["conditions"])
            self.assertEqual(pooled.EXPECTED_MEMBER_HISTOGRAM[source],
                             common.released_histogram(source))
            self.assertEqual(pooled.EXPECTED_ROWS[source], released["rows"])
            self.assertEqual(pooled.EXPECTED_TREATED[source], released["treated"])
            self.assertEqual(pooled.EXPECTED_VEHICLE[source], released["vehicle"])
            self.assertEqual(truth.EXPECTED_CONDITIONS[source], released["conditions"])
            self.assertEqual(truth.EXPECTED_REPLICATED[source],
                             common.replicated_conditions()[source])

    def test_released_shape_reconciles(self):
        """Rows split into treated+vehicle and the member counts add up."""
        for source in common.SOURCE_ORDER:
            released = common.RELEASED[source]
            self.assertEqual(released["treated"] + released["vehicle"], released["rows"])
            histogram = common.released_histogram(source)
            self.assertEqual(sum(int(count) * int(members)
                                 for members, count in histogram.items()),
                             released["treated"])


if __name__ == "__main__":
    unittest.main()
