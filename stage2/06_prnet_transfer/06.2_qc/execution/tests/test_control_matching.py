"""Control-matching tests: the join keys and the "never cross-line / never mix arms" rules."""

import unittest

import pandas as pd

import control_matching as cm


def _frame(rows):
    return pd.DataFrame(rows, index=[r["observation_id"] for r in rows])


def _tempenv():
    return _frame([
        # plate 1: one eligible DMSO, one QC-failed DMSO, a Media well we must NOT use
        dict(observation_id="temp_env::Plate1-B02", plate="Plate1", vehicle="DMSO",
             role="vehicle_control_dmso", qc_status="pass"),
        dict(observation_id="temp_env::Plate1-B03", plate="Plate1", vehicle="DMSO",
             role="vehicle_control_dmso", qc_status="fail"),
        dict(observation_id="temp_env::Plate1-B04", plate="Plate1", vehicle="Media",
             role="control_media", qc_status="pass"),
        dict(observation_id="temp_env::Plate1-B12", plate="Plate1", vehicle="Blank",
             role="technical_blank", qc_status="pass"),
        dict(observation_id="temp_env::Plate1-B05", plate="Plate1", vehicle="DMSO",
             role="treated_compound", qc_status="pass"),
        # plate 2: every DMSO well failed QC -> its treated well is unmatched
        dict(observation_id="temp_env::Plate2-B02", plate="Plate2", vehicle="DMSO",
             role="vehicle_control_dmso", qc_status="fail"),
        dict(observation_id="temp_env::Plate2-B05", plate="Plate2", vehicle="DMSO",
             role="treated_compound", qc_status="pass"),
    ])


class TestTempEnvMatching(unittest.TestCase):
    def test_keys_are_plate_and_vehicle(self):
        f = _tempenv()
        keys = cm.control_key(f, "temp_env")
        self.assertEqual(keys.iloc[0], "temp_env::Plate1::vehicle=DMSO")
        self.assertEqual(keys.iloc[2], "temp_env::Plate1::vehicle=Media")

    def test_treated_counts_only_qc_passing_dmso_on_its_plate(self):
        m = cm.match_controls(_tempenv(), "temp_env")
        row = m.loc["temp_env::Plate1-B05"]
        self.assertEqual(int(row["matched_control_n"]), 1)      # the Media well is not a control
        self.assertEqual(int(row["matched_control_n_total"]), 2)  # two DMSO wells before QC
        self.assertTrue(bool(row["matched_control_ok"]))
        self.assertEqual(row["matched_control_status"], "matched")

    def test_plate_with_no_passing_vehicle_is_unmatched(self):
        m = cm.match_controls(_tempenv(), "temp_env")
        row = m.loc["temp_env::Plate2-B05"]
        self.assertEqual(int(row["matched_control_n"]), 0)
        self.assertFalse(bool(row["matched_control_ok"]))
        self.assertEqual(row["matched_control_status"], "unmatched")

    def test_control_rows_carry_no_verdict(self):
        m = cm.match_controls(_tempenv(), "temp_env")
        self.assertTrue(pd.isna(m.loc["temp_env::Plate1-B02", "matched_control_n"]))
        self.assertEqual(m.loc["temp_env::Plate1-B02", "matched_control_status"], "is_control")

    def test_media_and_blank_are_not_arms_and_never_enter_the_pool(self):
        m = cm.match_controls(_tempenv(), "temp_env")
        for wid in ("temp_env::Plate1-B04", "temp_env::Plate1-B12"):
            self.assertTrue(pd.isna(m.loc[wid, "matched_control_n"]))
            self.assertEqual(m.loc[wid, "matched_control_status"], "not_an_arm")

    def test_availability_table_reports_unmatched_and_qc_loss(self):
        t = cm.availability_table(_tempenv(), "temp_env")
        p1 = t[(t.plate == "Plate1") & (t.vehicle == "DMSO")].iloc[0]
        self.assertEqual(int(p1["control_observations_total"]), 2)
        self.assertEqual(int(p1["control_observations_qc_pass"]), 1)
        self.assertEqual(int(p1["treated_observations_qc_pass"]), 1)
        self.assertEqual(int(p1["treated_kept_unmatched"]), 0)
        p2 = t[(t.plate == "Plate2") & (t.vehicle == "DMSO")].iloc[0]
        self.assertEqual(int(p2["treated_kept_unmatched"]), 1)


def _tahoe():
    return _frame([
        # plate1 / line A: one passing control, one failing control
        dict(observation_id="t1", plate="plate1", cell_line_id="CVCL_A", role="vehicle_control_dmso_tf",
             qc_status="pass"),
        dict(observation_id="t2", plate="plate1", cell_line_id="CVCL_A", role="vehicle_control_dmso_tf",
             qc_status="fail"),
        dict(observation_id="t3", plate="plate1", cell_line_id="CVCL_A", role="treated_compound",
             qc_status="pass"),
        # plate1 / line B has a control, but it must never serve line A
        dict(observation_id="t4", plate="plate1", cell_line_id="CVCL_B", role="vehicle_control_dmso_tf",
             qc_status="pass"),
        dict(observation_id="t5", plate="plate1", cell_line_id="CVCL_B", role="treated_compound",
             qc_status="pass"),
        # plate2 / line A: both controls failed
        dict(observation_id="t6", plate="plate2", cell_line_id="CVCL_A", role="vehicle_control_dmso_tf",
             qc_status="fail"),
        dict(observation_id="t7", plate="plate2", cell_line_id="CVCL_A", role="treated_compound",
             qc_status="pass"),
    ])


class TestTahoeMatching(unittest.TestCase):
    def test_key_is_plate_and_cell_line(self):
        keys = cm.control_key(_tahoe(), "tahoe")
        self.assertEqual(keys.loc["t3"], "tahoe::plate1::cell_line_id=CVCL_A")
        self.assertNotEqual(keys.loc["t3"], keys.loc["t5"])

    def test_never_crosses_cell_line(self):
        m = cm.match_controls(_tahoe(), "tahoe")
        self.assertEqual(int(m.loc["t3", "matched_control_n"]), 1)  # CVCL_B's control does not count
        self.assertTrue(bool(m.loc["t3", "matched_control_ok"]))

    def test_group_with_all_controls_failing_is_unmatched(self):
        m = cm.match_controls(_tahoe(), "tahoe")
        self.assertEqual(int(m.loc["t7", "matched_control_n"]), 0)
        self.assertEqual(m.loc["t7", "matched_control_status"], "unmatched")

    def test_availability_table_counts_unmatched_line_groups(self):
        t = cm.availability_table(_tahoe(), "tahoe")
        self.assertEqual(len(t), 3)
        self.assertEqual(int(t.treated_kept_unmatched.sum()), 1)


class TestKeyCoverage(unittest.TestCase):
    def test_every_source_has_both_a_key_and_a_control_role(self):
        for source in ("temp_env", "tahoe"):
            self.assertIn(source, cm.CONTROL_KEYS)
            self.assertIn(source, cm.CONTROL_ROLES)


if __name__ == "__main__":
    unittest.main()
