"""Pinned crosswalk and class-rule tests for the cardiotoxicity sources (#189).

Every fixture is hand-built and small: the tests assert the *rule* (annotation
resolution, the ``(lookup_name, dataset_inchikey)`` key, the ATC level-4
partition, TOP2 corroboration, the >=3 small-molecule rule and the recorded
exclusions), not a re-derivation from the implementation's internals.
"""
import gzip
import io
import json
import os
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import cardiotox, disease_deg, release, r_env

TOP2 = "Topoisomerase II Inhibitors"
TKI = "Tyrosine Kinase Inhibitors"


def annotation_rows():
    return [
        # lookup_name, dataset_inchikey, parent_inchikey, atc, moa, mesh_classes
        ("Doxorubicin", "D1", "D1", "L01DB01 Anthracyclines and related substances",
         "", "Antibiotics, Antineoplastic; " + TOP2),
        # Name-only duplicate with a different class: the dataset row must win.
        ("Doxorubicin", "", "OTHER", "C01BD01 Antiarrhythmics, class III", "",
         "Anti-Arrhythmia Agents"),
        ("Daunorubicin", "D2", "D2", "L01DB02 Anthracyclines and related substances",
         "DNA topoisomerase II alpha inhibitor", TOP2),
        ("Epirubicin", "D4", "D4", "L01DB03 Anthracyclines and related substances",
         "", TOP2),
        ("Idarubicin", "D5", "D5", "L01DB06 Anthracyclines and related substances",
         "", TOP2),
        ("Mitoxantrone", "MTX", "MTX", "L01DB07 Anthracyclines and related substances",
         "", TOP2),
        # A TOP2-subgroup member without TOP2 mesh evidence: recorded excluded.
        ("Valrubicin", "VA", "VA", "L01DB09 Anthracyclines and related substances",
         "", "Antibiotics, Antineoplastic"),
        ("Afatinib", "AF", "AF", "L01EB03 Epidermal growth factor receptor kinase",
         "", TKI),
        ("Erlotinib", "E1", "E1", "L01EB02 Epidermal growth factor receptor kinase",
         "", TKI),
        ("Gefitinib", "G1", "G1", "L01EB01 Epidermal growth factor receptor kinase",
         "", TKI),
        ("Dabrafenib", "DAB", "DAB", "L01EC02 B-Raf serine-threonine kinase",
         "", "Protein Kinase Inhibitors"),
        ("Vemurafenib", "VEM", "VEM", "L01EC01 B-Raf serine-threonine kinase",
         "", "Protein Kinase Inhibitors"),
        # Multi-ATC code whose mesh evidence matches neither subgroup.
        ("Cyclosporine", "CY", "CY",
         "S01XA18 Other ophthalmologicals; L04AD01 Other immunosuppressants", "",
         "Antifungal Agents"),
        # Level-4-only ATC string (no level-5 suffix).
        ("Ruxolitinib", "RUX", "RUX", "D11AH Agents for dermatitis", "",
         "Janus Kinase Inhibitors"),
        # The DToxS AMI analogue: the dataset structure has no class evidence and
        # must not inherit the real amiodarone's.
        ("Amiodarone", "A1", "A1", "", "", ""),
        ("Amiodarone", "", "REAL", "C01BD01 Antiarrhythmics, class III", "",
         "Anti-Arrhythmia Agents"),
    ]


def write_annotation(root):
    path = root / "cm_compound_annotation.csv"
    pd.DataFrame(annotation_rows(), columns=[
        "lookup_name", "dataset_inchikey", "parent_inchikey", "atc", "moa",
        "mesh_classes"]).assign(
        inchikey="", smiles="", pubchem_cid="", chembl_id="", verification="verified",
        notes="").to_csv(path, index=False)
    return path


def write_series(root, name, pairs):
    """A minimal GEO series matrix with treatment/drug-name characteristic rows."""
    path = root / name
    with gzip.open(path, "wt") as handle:
        handle.write("!Sample_title\n")
        conditions = "\t".join(f'"{code} MSN01-CM Rep 1"' for code in pairs)
        handle.write("!Sample_title\t" + conditions + "\n")
        handle.write("!Sample_characteristics_ch1\t"
                     + "\t".join(f'"treatment condition: {code}"' for code in pairs)
                     + "\n")
        handle.write("!Sample_characteristics_ch1\t"
                     + "\t".join(f'"drug name: {pairs[code]}"'
                                 for code in pairs) + "\n")
    return path


def write_dtoxs_params(root, codes):
    path = root / "GSE217421_Conv-RNAseq-Params.tar.gz"
    payload = "Group\tDrug\tType\n" + "".join(
        f"{code}\t{code}\tOffender\n" for code in codes)
    with tarfile.open(path, "w:gz") as archive:
        data = payload.encode()
        info = tarfile.TarInfo("Conv-RNAseq-Drug-Groups.tsv")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    return path


def write_gse244740_series(root):
    wells = ["S4881Nr8_A1", "S4881Nr8_A2", "S4881Nr8_A3", "S4881Nr8_A4"]
    compounds = ["Doxorubicin", "Idarubicin", "Mitoxantrone", "DMSO"]
    path = root / "GSE244740_series_matrix.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("!Sample_description\t" + "\t".join(f'"{well}"' for well in wells) + "\n")
        handle.write("!Sample_characteristics_ch1\t"
                     + "\t".join(f'"treatment: {compound}, Con8"'
                                 for compound in compounds) + "\n")
        handle.write("!Sample_characteristics_ch1\t"
                     + "\t".join('"time: 24h"' for _ in wells) + "\n")
    return path


def write_gse243674_json(root):
    path = root / "GSE243674_deg.json"
    path.write_text(json.dumps({"contrasts": [
        {"contrast": "dnr_3h_vs_vehicle", "drug_token": "DNR", "drug": "daunorubicin",
         "time_h": 3, "n_case": 6, "n_control": 6},
        {"contrast": "dox_3h_vs_vehicle", "drug_token": "DOX", "drug": "doxorubicin",
         "time_h": 3, "n_case": 6, "n_control": 6},
        {"contrast": "epi_3h_vs_vehicle", "drug_token": "EPI", "drug": "epirubicin",
         "time_h": 3, "n_case": 6, "n_control": 6},
        {"contrast": "mtx_3h_vs_vehicle", "drug_token": "MTX", "drug": "mitoxantrone",
         "time_h": 3, "n_case": 6, "n_control": 6},
        {"contrast": "trz_3h_vs_vehicle", "drug_token": "TRZ", "drug": "trastuzumab",
         "time_h": 3, "n_case": 6, "n_control": 6},
    ]}))
    return path


DTOXS_CODES = {"DAU": "Daunorubicin", "DOX": "Doxorubicin", "EPI": "Epirubicin",
               "IDA": "Idarubicin", "VAL": "Valrubicin", "AFA": "Afatinib",
               "ERL": "Erlotinib", "GEF": "Gefitinib", "DAB": "Dabrafenib",
               "VEM": "Vemurafenib", "CYC": "Cyclosporine", "RUX": "Ruxolitinib",
               "AMI": "Amiodarone", "TRS": "Trastuzumab"}


class Fixture:
    def __init__(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "pert_data/DToxS-iPSC/source").mkdir(parents=True)
        (self.root / "pert_data/GSE244740/source").mkdir(parents=True)
        (self.root / "disease_data/GSE243674/local/deg").mkdir(parents=True)
        (self.root / "pert_data/Annotation/local").mkdir(parents=True)
        self.annotation = write_annotation(self.root / "pert_data/Annotation/local")
        write_series(self.root / "pert_data/DToxS-iPSC/source",
                     "GSE217421_series_matrix.txt.gz", DTOXS_CODES)
        write_dtoxs_params(self.root / "pert_data/DToxS-iPSC/source",
                           sorted(DTOXS_CODES))
        write_gse244740_series(self.root / "pert_data/GSE244740/source")
        write_gse243674_json(self.root / "disease_data/GSE243674/local/deg")

    def close(self):
        self.temporary.cleanup()


class CrosswalkTest(unittest.TestCase):
    def setUp(self):
        self.fixture = Fixture()
        self.root = self.fixture.root
        self.annotation = cardiotox.read_annotation(self.fixture.annotation)

    def tearDown(self):
        self.fixture.close()

    def crosswalk(self):
        return cardiotox.build_crosswalk(self.root, self.annotation)

    def test_source_records_cover_every_source_code(self):
        records = cardiotox.source_records(self.root)
        by_source = {}
        for record in records:
            by_source.setdefault(record["source"], set()).add(record["source_id"])
        self.assertEqual(by_source["dtoxs_ipsc"], set(DTOXS_CODES))
        self.assertEqual(by_source["gse243674"], {"DNR", "DOX", "EPI", "MTX", "TRZ"})
        self.assertEqual(by_source["gse244740"],
                         {"Doxorubicin", "Idarubicin", "Mitoxantrone"})

    def test_dataset_inchikey_wins_over_the_name_only_row(self):
        crosswalk = self.crosswalk()
        dox = crosswalk[(crosswalk.source == "dtoxs_ipsc")
                        & (crosswalk.source_id == "DOX")].iloc[0]
        self.assertEqual(dox.dataset_inchikey, "D1")
        self.assertEqual(dox.class_id, "atc_l4:L01DB")

    def test_amiodarone_analogue_is_not_relabelled(self):
        crosswalk = self.crosswalk()
        ami = crosswalk[(crosswalk.source == "dtoxs_ipsc")
                        & (crosswalk.source_id == "AMI")].iloc[0]
        self.assertEqual(ami.dataset_inchikey, "A1")
        self.assertEqual(ami.lookup_name, "Amiodarone")
        self.assertEqual(ami.status, "unresolved_no_class")
        self.assertEqual(ami.class_id, "")

    def test_missing_annotation_is_recorded_not_guessed(self):
        crosswalk = self.crosswalk()
        trs = crosswalk[crosswalk.source_id == "TRS"].iloc[0]
        self.assertEqual(trs.status, "unresolved_no_annotation")
        self.assertFalse(bool(trs.small_molecule))
        self.assertIn("no annotation row", trs.exclusion_reason)

    def test_multi_atc_without_mesh_corroboration_is_ambiguous(self):
        crosswalk = self.crosswalk()
        cyc = crosswalk[(crosswalk.source == "dtoxs_ipsc")
                        & (crosswalk.source_id == "CYC")].iloc[0]
        self.assertEqual(cyc.status, "ambiguous_multi_class")
        classification = cardiotox.classify_row(
            {"atc": "S01XA18 Other ophthalmologicals; L04AD01 Other immunosuppressants",
             "mesh_classes": "Antifungal Agents"})
        self.assertEqual(classification["status"], "ambiguous_multi_class")
        self.assertEqual(classification["atc_candidates"],
                         ["atc_l4:L04AD", "atc_l4:S01XA"])

    def test_level4_only_atc_string_is_a_class(self):
        entries = cardiotox.atc_entries("D11AH Agents for dermatitis")
        self.assertEqual(entries, [("D11AH", "Agents for dermatitis")])
        self.assertEqual(cardiotox.atc_level4("D11AH"), "D11AH")
        crosswalk = self.crosswalk()
        rux = crosswalk[(crosswalk.source == "dtoxs_ipsc")
                        & (crosswalk.source_id == "RUX")].iloc[0]
        self.assertEqual(rux.class_id, "atc_l4:D11AH")
        self.assertTrue(bool(rux.small_molecule))

    def test_partition_every_assigned_code_has_one_class(self):
        crosswalk = self.crosswalk()
        assigned = crosswalk[crosswalk.status == "assigned"]
        self.assertTrue(assigned.class_id.ne("").all())
        self.assertTrue(assigned.class_id.str.startswith("atc_l4:").all())

    def test_crosswalk_is_deterministic(self):
        pd.testing.assert_frame_equal(self.crosswalk(), self.crosswalk())


class ClassRuleTest(unittest.TestCase):
    def setUp(self):
        self.fixture = Fixture()
        self.root = self.fixture.root
        annotation = cardiotox.read_annotation(self.fixture.annotation)
        self.crosswalk = cardiotox.build_crosswalk(self.root, annotation)
        self.classes, self.members, self.excluded = cardiotox.select_classes(
            self.crosswalk)

    def tearDown(self):
        self.fixture.close()

    def _class(self, source, class_id):
        return self.classes[(self.classes.source == source)
                            & (self.classes.class_id == class_id)].iloc[0]

    def test_three_small_molecule_rule(self):
        self.assertEqual(self._class("dtoxs_ipsc", "atc_l4:L01EB").n_small_molecule_drugs, 3)
        self.assertEqual(self._class("dtoxs_ipsc", "atc_l4:L01EB").decision, "entered")
        below = self._class("dtoxs_ipsc", "atc_l4:L01EC")
        self.assertEqual(below.n_small_molecule_drugs, 2)
        self.assertEqual(below.decision, "not_entered")
        self.assertEqual(below.reason, "class_below_three")

    def test_top2_class_enters_each_source_and_other_only_from_dtoxs(self):
        entered = self.classes[self.classes.decision == "entered"]
        self.assertEqual(
            set(entered[entered.category == "cardiotox_TOP2"].source),
            {"dtoxs_ipsc", "gse243674", "gse244740"})
        self.assertEqual(
            set(entered[entered.category == "cardiotox_other"].class_id),
            {"atc_l4:L01EB"})

    def test_top2_mesh_missing_member_is_recorded_excluded(self):
        val = self.excluded[(self.excluded.source == "dtoxs_ipsc")
                            & (self.excluded.source_id == "VAL")]
        self.assertEqual(len(val), 1)
        self.assertEqual(val.iloc[0].reason, "top2_mesh_missing")
        members = self.members[(self.members.source == "dtoxs_ipsc")
                               & (self.members.class_id == "atc_l4:L01DB")]
        self.assertNotIn("VAL", set(members.source_id))
        self.assertEqual(self._class("dtoxs_ipsc", "atc_l4:L01DB").n_small_molecule_drugs, 4)

    def test_antibody_and_unresolved_codes_are_excluded(self):
        excluded = set(self.excluded.source_id)
        self.assertIn("TRS", excluded)
        for name in ("Trastuzumab",):
            self.assertIn(name, set(self.crosswalk[self.crosswalk.status
                                                   == "unresolved_no_annotation"]
                                    .drug_name))
        self.assertEqual(
            self.crosswalk[(self.crosswalk.source_id == "TRS")
                           & (self.crosswalk.small_molecule.astype(bool))].shape[0], 0)

    def test_every_non_member_code_is_recorded_in_exclusions(self):
        member_keys = {(row.source, row.source_id)
                       for row in self.members.itertuples(index=False)}
        excluded_keys = {(row.source, row.source_id)
                         for row in self.excluded.itertuples(index=False)}
        for row in self.crosswalk.itertuples(index=False):
            if (row.source, row.source_id) not in member_keys:
                self.assertIn((row.source, row.source_id), excluded_keys,
                              f"{row.source}:{row.source_id} not recorded as excluded")

    def test_class_counts_are_per_source_and_across_sources(self):
        row = self._class("dtoxs_ipsc", "atc_l4:L01DB")
        self.assertEqual(row.n_small_molecule_drugs, 4)
        # Doxorubicin/daunorubicin-like names repeat across sources, so the
        # across-source count is a distinct-name count over resolved members.
        self.assertGreaterEqual(int(row.n_small_molecule_across_sources), 4)


class EvidenceProductTest(unittest.TestCase):
    def setUp(self):
        self.fixture = Fixture()
        self.root = self.fixture.root
        annotation = cardiotox.read_annotation(self.fixture.annotation)
        self.crosswalk = cardiotox.build_crosswalk(self.root, annotation)
        self.classes, self.members, self.excluded = cardiotox.select_classes(
            self.crosswalk)

    def tearDown(self):
        self.fixture.close()

    def test_evidence_product_is_deterministic_and_relative(self):
        first = Path(tempfile.mkdtemp())
        second = Path(tempfile.mkdtemp())
        try:
            manifests = []
            for target in (first, second):
                manifests.append(cardiotox.write_evidence_product(
                    target / "cardiotox", self.root, self.crosswalk, self.classes,
                    self.members, self.excluded, cardiotox.rejected_table(),
                    {}, {"python": "3.13"}, "0" * 64))
            for name in ("crosswalk.csv", "classes.csv", "class_members.csv",
                         "excluded.csv", "eligibility.csv", "rejected.csv",
                         "manifest.json"):
                self.assertEqual(release.sha256(first / "cardiotox" / name),
                                 release.sha256(second / "cardiotox" / name), name)
            release.check_relative_paths(manifests[0])
            release.check_manifest(first / "cardiotox", manifests[0])
            self.assertIn("crosswalk.csv", manifests[0]["outputs"])
        finally:
            import shutil
            shutil.rmtree(first, ignore_errors=True)
            shutil.rmtree(second, ignore_errors=True)

    def test_eligibility_records_every_source(self):
        eligibility = cardiotox.eligibility_table(
            self.root, self.crosswalk, self.classes, self.members, self.excluded)
        self.assertEqual(set(eligibility.source), set(cardiotox.SOURCES))
        dtoxs = eligibility[eligibility.source == "dtoxs_ipsc"].iloc[0]
        self.assertEqual(dtoxs.evidence_type, "bulk")
        self.assertEqual(int(dtoxs.tier), 2)
        self.assertIn("GSE217423", dtoxs.notes)
        categories = set(self.classes[(self.classes.source == "dtoxs_ipsc")
                                      & (self.classes.decision == "entered")]
                         .category)
        self.assertEqual(categories, {"cardiotox_TOP2", "cardiotox_other"})

    def test_rejected_sources_record_the_platform_and_single_drug_series(self):
        rejected = cardiotox.rejected_table()
        text = " ".join(rejected.candidate + " " + rejected.reason)
        self.assertIn("GSE217423", text)
        self.assertIn("GSE242692", text)
        self.assertIn("GSE313746", text)


class RResolutionTests(unittest.TestCase):
    """The R launchers resolve configuration, never a host-absolute literal (#194).

    No R is required: the overrides point at throwaway shell scripts, and the
    missing-interpreter case is exercised with an empty ``PATH``.
    """

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def fake_rscript(self, body="exit 0\n", name="Rscript"):
        path = self.root / name
        path.write_text("#!/bin/sh\n" + body)
        path.chmod(0o755)
        return path

    def test_arc_rscript_override_wins_for_the_cardiotox_launcher(self):
        marker = self.root / "ran"
        fake = self.fake_rscript(body=f": > {marker}\nexit 0\n")
        with mock.patch.dict(os.environ, {"ARC_RSCRIPT": str(fake)}, clear=True):
            self.assertEqual(r_env.rscript(), str(fake))
            cardiotox.run_r("step.R", [])
        self.assertTrue(marker.is_file())

    def test_arc_r_home_override_and_rscript_derivation_for_the_deg_worker(self):
        home = self.root / "lib" / "R"
        (home / "bin" / "exec").mkdir(parents=True)
        (home / "bin" / "exec" / "R").write_text("")
        with mock.patch.dict(os.environ, {"ARC_R_HOME": str(home)}, clear=True):
            executable, resolved = disease_deg.r_entry_point()
        self.assertEqual((executable, resolved), (home / "bin" / "exec" / "R", home))
        printed = self.fake_rscript(body=f"echo {home}\nexit 0\n", name="print_home")
        bare_env = mock.patch.object(sys, "prefix", str(self.root / "no-env-r"))
        with mock.patch.dict(os.environ, {"ARC_RSCRIPT": str(printed)}, clear=True), bare_env:
            self.assertEqual(r_env.r_home(), home)
            self.assertEqual(disease_deg.r_entry_point(),
                             (home / "bin" / "exec" / "R", home))

    def test_missing_rscript_raises_a_clear_error(self):
        empty = self.root / "empty"
        empty.mkdir()
        bare_env = mock.patch.object(sys, "prefix", str(self.root / "no-env-r"))
        with mock.patch.dict(os.environ, {"PATH": str(empty)}, clear=True), bare_env:
            with self.assertRaisesRegex(r_env.REnvironmentError, "Rscript"):
                r_env.rscript()
            with self.assertRaisesRegex(r_env.REnvironmentError, "Rscript"):
                disease_deg.r_entry_point()
            with self.assertRaisesRegex(r_env.REnvironmentError, "Rscript"):
                cardiotox.run_r("step.R", [])

    def test_the_running_env_r_wins_over_a_bare_path(self):
        bundled = Path(sys.prefix) / "lib" / "R"
        if not (bundled / "bin" / "exec" / "R").is_file():
            self.skipTest("the running conda env ships no R")
        empty = self.root / "empty"
        empty.mkdir()
        with mock.patch.dict(os.environ, {"PATH": str(empty)}, clear=True):
            self.assertEqual(r_env.r_home(), bundled)


if __name__ == "__main__":
    unittest.main()
