"""Ticket 03 item 4: the tracked per-condition discarded-member evidence table.

The published member-choice payloads carry a ``members.csv`` companion, but that file is
gitignored with the rest of the payload, so the number of discarded members per condition
was only derivable. ``log/ticket-03-discarded-members.csv`` makes it a reviewable tracked
product. This test does two things:

* always: reconcile the tracked table with the released-population constants (row counts,
  the member-count histogram and the discarded totals);
* with the published payloads present: recompute the same counts from the *merge* member
  table, an independently built product, and regenerate the tracked table byte-for-byte.
"""
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import common  # noqa: E402
import evidence  # noqa: E402
import tsg_view_fixtures as fixtures  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
EVIDENCE = HERE.parents[1] / "log/ticket-03-discarded-members.csv"
VIEWS = fixtures.payload_root(REPO) / "stage2/tsg_view/output/views"
#: ``solo`` is published for both sources, ``solo-deepest`` (the sensitivity arm) for tenv.
EXPECTED_ROWS = {"draw": {"tenv": 1377, "tahoe_bulk": 52777},
                 "deepest": {"tenv": 1377}}
NEEDS_PAYLOADS = unittest.skipUnless(
    (VIEWS / "merge/tenv.members.csv").is_file()
    and (VIEWS / "solo/tenv.members.csv").is_file(),
    fixtures.dependency_message("the published view payloads", VIEWS))


def read_evidence():
    return pd.read_csv(EVIDENCE, dtype={"profile_id": str, "chosen_member_id": str,
                                        "selection_rule": str, "source": str})


def histogram(frame):
    return {str(key): int(value)
            for key, value in frame["member_count"].value_counts().sort_index().items()}


class TrackedEvidenceTest(unittest.TestCase):
    """The committed table, checked without touching a payload."""

    def test_columns_and_one_row_per_condition(self):
        table = read_evidence()
        self.assertEqual(tuple(table.columns), evidence.COLUMNS)
        self.assertEqual(len(table),
                         sum(sum(block.values()) for block in EXPECTED_ROWS.values()))
        duplicate = table.duplicated(["source", "selection_rule", "profile_id"])
        self.assertFalse(duplicate.any(), "a condition is listed twice")
        self.assertTrue((table["profile_id"] == table["chosen_member_id"]).all())

    def test_counts_reconcile_with_the_released_population(self):
        table = read_evidence()
        for rule, expected in EXPECTED_ROWS.items():
            for source, rows in expected.items():
                block = table[(table["selection_rule"] == rule)
                              & (table["source"] == source)]
                self.assertEqual(len(block), rows)
                self.assertTrue((block["discarded_count"]
                                 == block["member_count"] - 1).all())
                released = {key: value for key, value in
                            common.released_histogram(source).items() if value}
                self.assertEqual(histogram(block), released)
                self.assertEqual(int(block["discarded_count"].sum()),
                                 sum((int(members) - 1) * count
                                     for members, count in released.items()))
                if source == "tenv":
                    self.assertEqual(histogram(block), {"3": 1320, "4": 48, "6": 9})
                else:
                    self.assertEqual(histogram(block)["1"], 43036)
                    self.assertEqual(histogram(block)["16"], 50)

    def test_the_two_tenv_arms_share_the_population_but_not_every_choice(self):
        table = read_evidence()
        draw = table[(table["selection_rule"] == "draw") & (table["source"] == "tenv")]
        deepest = table[(table["selection_rule"] == "deepest")
                        & (table["source"] == "tenv")]
        self.assertEqual(histogram(draw), histogram(deepest))
        self.assertEqual(sorted(draw["member_count"]), sorted(deepest["member_count"]))
        # The rules differ: on a multi-member condition they may keep another member.
        self.assertNotEqual(sorted(draw["profile_id"]), sorted(deepest["profile_id"]))


class RegenerationTest(unittest.TestCase):
    def test_a_missing_member_table_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "member table is missing"):
                evidence.discarded_members(tmp)

    @NEEDS_PAYLOADS
    def test_regeneration_is_byte_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, frame = evidence.write_evidence(VIEWS, Path(tmp) / "table.csv")
            self.assertEqual(path.read_bytes(), EVIDENCE.read_bytes())
            self.assertEqual(len(frame), len(read_evidence()))


class IndependentRecomputeTest(unittest.TestCase):
    """Recompute every count from the merge member table, a different product."""

    @NEEDS_PAYLOADS
    def test_counts_match_the_merge_member_table(self):
        table = read_evidence()
        for source in common.SOURCE_ORDER:
            merge = pd.read_csv(VIEWS / "merge" / f"{source}.members.csv", dtype=str)
            # The merge member table lists every member of every condition under that
            # condition's pooled identity, so it yields an independent member count per
            # profile without reading any member-choice product.
            sizes = merge.groupby("profile_id", sort=False).size()
            member_count = {}
            for identity, block in merge.groupby("profile_id", sort=False):
                for member in block["member_profile_id"]:
                    member_count[member] = len(block)
            for rule in table[table["source"] == source]["selection_rule"].unique():
                observed = table[(table["source"] == source)
                                 & (table["selection_rule"] == rule)]
                self.assertTrue(set(observed["profile_id"]) <= set(member_count),
                                f"{source}: a chosen member is not in the merge member table")
                for row in observed.itertuples():
                    self.assertEqual(int(row.member_count), member_count[row.profile_id])
                # Both arms discard the same number per condition: the population is the
                # merge table's, only the kept member differs.
                self.assertEqual(int(observed["discarded_count"].sum()),
                                 int((sizes - 1).sum()))

    @NEEDS_PAYLOADS
    def test_a_singleton_condition_discards_nothing(self):
        table = read_evidence()
        singleton = table[table["member_count"] == 1]
        self.assertEqual(len(singleton), 43036)
        self.assertTrue((singleton["discarded_count"] == 0).all())
        self.assertTrue((singleton["selection_rule"] == "draw").all())
        self.assertTrue((singleton["source"] == "tahoe_bulk").all())


if __name__ == "__main__":
    unittest.main()
