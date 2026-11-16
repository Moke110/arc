"""Accepted ADR 0072 five-category product contracts at the public file boundary.

Non-TOX records retain bare input digests and item-local derived inputs; TOX
pins measured source inputs directly as {bytes, sha256} records.
"""
import hashlib
import json
import re
import unittest
from pathlib import Path

REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
DISEASE_DEG = REPO / "stage2" / "01_bases" / "output" / "disease-deg"
EXPECTED_CATEGORIES = (
    "af",
    "icm",
    "nhf",
    "genetic_cardiomyopathy",
    "tox",
)
NON_TOX = {"af", "icm", "nhf", "genetic_cardiomyopathy"}
SCHEMA = "stage2-disease-deg-catalog-v1"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
RECORD_FIELDS = {
    "dataset_id", "category", "evidence_type", "tier", "deg", "rank",
    "method", "inputs", "n_measured", "n_significant",
    # ADR 0068 / #194: every record states whether its measured effect scores
    # carry an ordering, so the fusion never has to guess it.
    "informative", "informativeness",
}
LEGACY_CATALOG_KEYS = {"test_method", "dataset_records"}
LEGACY_RECORD_KEYS = {"test_method", "dataset_records"}


def catalog_digest(catalog):
    """ADR 0062/#189 canonical digest: content without the self-reference."""
    body = dict(catalog)
    body.pop("catalog_sha256", None)
    return hashlib.sha256(
        (json.dumps(body, indent=2, sort_keys=True) + "\n").encode()).hexdigest()


class CatalogShapeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalogs = {}
        for path in sorted(DISEASE_DEG.glob("*/catalog.json")):
            cls.catalogs[path.parent.name] = json.loads(path.read_text())

    def test_all_categories_have_catalogs_and_nothing_else(self):
        self.assertEqual(set(self.catalogs), set(EXPECTED_CATEGORIES))
        self.assertEqual(len(self.catalogs), len(EXPECTED_CATEGORIES))


    def test_catalog_schema_and_records(self):
        for category, catalog in self.catalogs.items():
            with self.subTest(category=category):
                self.assertEqual(catalog.get("schema"), SCHEMA)
                self.assertFalse(LEGACY_CATALOG_KEYS & set(catalog),
                                 f"{category}: retired catalog keys present")
                datasets = catalog.get("datasets")
                self.assertIsInstance(datasets, list)
                self.assertGreater(len(datasets), 0)
                for record in datasets:
                    self.assertIsInstance(record, dict)
                    self.assertFalse(LEGACY_RECORD_KEYS & set(record),
                                     f"{category}: retired record keys present")
                    self.assertTrue(RECORD_FIELDS <= set(record),
                                    f"{category}: record fields missing: "
                                    f"{sorted(RECORD_FIELDS - set(record))}")
                    self.assertIsInstance(record["method"], dict)
                    self.assertTrue(record["method"],
                                    "a dataset record must carry its method")
                    self.assertTrue(str(record["dataset_id"]))

    def test_informativeness_is_explicit_and_reasoned(self):
        # ADR 0068 / ADR 0071: a dataset with no v2 DEG (or all-tied effects) is
        # published with informative: false and a recorded reason; the field is a
        # boolean, never inferred downstream. A zero-DEG dataset may still have
        # many distinct effect scores -- its ordering evidence is absent because
        # nothing is a DEG, not because the scores are tied -- so only the
        # all-tied case is required to have n_distinct_effect == 1.
        for category, catalog in self.catalogs.items():
            for record in catalog["datasets"]:
                with self.subTest(category=category, dataset=record["dataset_id"]):
                    self.assertIsInstance(record["informative"], bool)
                    details = record["informativeness"]
                    if category != "tox":
                        self.assertIsInstance(details["rule"], str)
                    self.assertGreaterEqual(int(details["n_distinct_effect"]), 1)
                    if record["informative"]:
                        self.assertIsNone(details["reason"])
                        self.assertGreater(details["n_distinct_effect"], 1)
                    else:
                        self.assertIsInstance(details["reason"], str)
                        self.assertTrue(details["reason"])
                        if int(record.get("n_deg", 0)) == 0:
                            if category != "tox":
                                self.assertIn("none of the", details["reason"])
                        else:
                            self.assertEqual(int(details["n_distinct_effect"]), 1)

    def test_catalog_sha256_is_self_consistent_and_in_the_manifest(self):
        for category, catalog in self.catalogs.items():
            with self.subTest(category=category):
                declared = catalog.get("catalog_sha256")
                self.assertTrue(isinstance(declared, str) and HEX64.match(declared),
                                f"{category}: catalog_sha256 is not 64-hex")
                self.assertEqual(declared, catalog_digest(catalog),
                                 f"{category}: catalog_sha256 does not match content")
                manifest = json.loads(
                    (DISEASE_DEG / category / "manifest.json").read_text())
                self.assertEqual(manifest.get("catalog_sha256"), declared)

    def test_case_and_control_are_positive_or_have_a_declared_basis(self):
        for category, catalog in self.catalogs.items():
            for record in catalog["datasets"]:
                with self.subTest(category=category, dataset=record["dataset_id"]):
                    for side in ("n_case", "n_control"):
                        value = record.get(side)
                        basis = record.get(f"{side}_basis")
                        if isinstance(value, int) and value > 0:
                            continue
                        self.assertTrue(
                            isinstance(basis, str) and basis.strip(),
                            f"{record['dataset_id']}: {side}={value!r} without an "
                            f"explicit {side}_basis")

    def test_inputs_shape_variance_is_explicit(self):
        """The non-TOX bare-digest and measured TOX structured shapes, pinned."""
        for category, catalog in self.catalogs.items():
            expected_shape = "bare" if category in NON_TOX else "dict"
            for record in catalog["datasets"]:
                with self.subTest(category=category, dataset=record["dataset_id"]):
                    inputs = record.get("inputs")
                    self.assertIsInstance(inputs, dict)
                    self.assertGreater(len(inputs), 0)
                    for path, value in inputs.items():
                        self.assertFalse(Path(path).is_absolute(),
                                         f"absolute input path: {path}")
                        if expected_shape == "bare":
                            self.assertTrue(
                                isinstance(value, str) and HEX64.match(value),
                                f"{category}: expected a bare sha256 digest for "
                                f"{path}, got {type(value).__name__}")
                        else:
                            self.assertIsInstance(value, dict)
                            self.assertEqual(set(value), {"bytes", "sha256"})
                            self.assertIsInstance(value["bytes"], int)
                            self.assertTrue(HEX64.match(value["sha256"]))

    def test_derived_inputs_are_digest_pinned_where_declared(self):
        for category, catalog in self.catalogs.items():
            for record in catalog["datasets"]:
                with self.subTest(category=category, dataset=record["dataset_id"]):
                    derived = record.get("derived_inputs")
                    if category in NON_TOX:
                        self.assertIsInstance(derived, list)
                        self.assertGreater(len(derived), 0)
                        for entry in derived:
                            self.assertEqual(set(entry), {"relative", "sha256"})
                            self.assertFalse(Path(entry["relative"]).is_absolute())
                            self.assertTrue(HEX64.match(entry["sha256"]),
                                            f"{record['dataset_id']}: bad digest")
                    else:
                        # TOX pins original measurements directly in `inputs`.
                        self.assertIsNone(derived)

    def test_axis_declarations_resolve_to_the_same_frozen_axis(self):
        """Non-TOX pins the published copy, TOX the source table; both agree."""
        published = REPO / "stage2" / "01_bases" / "output" / "gene-identity" / "genes.csv"
        source = (REPO / "stage2" / "01_bases" / "output" / "gene-identity"
                  / "tables" / "genes.csv")
        published_digest = hashlib.sha256(published.read_bytes()).hexdigest()
        source_digest = hashlib.sha256(source.read_bytes()).hexdigest()
        self.assertNotEqual(published_digest, source_digest,
                            "the two declared axis files are the documented variance")
        for category, catalog in self.catalogs.items():
            with self.subTest(category=category):
                axis = catalog.get("axis")
                self.assertIsInstance(axis, dict)
                self.assertEqual(axis.get("size"), 18164)
                path = REPO / axis["path"]
                self.assertTrue(path.is_file(), f"{category}: axis path missing")
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                                 axis["sha256"])
                expected = (published_digest if category in NON_TOX
                            else source_digest)
                self.assertEqual(axis["sha256"], expected,
                                 f"{category}: unrecorded axis-declaration drift")

    def test_both_axis_files_are_the_same_18164_gene_order(self):
        published = (REPO / "stage2" / "01_bases" / "output" / "gene-identity"
                     / "genes.csv").read_text().splitlines()
        source = (REPO / "stage2" / "01_bases" / "output" / "gene-identity"
                  / "tables" / "genes.csv").read_text().splitlines()
        published_symbols = [line.split(",")[1] for line in published[1:]]
        source_symbols = [line.split(",")[0] for line in source[1:]]
        self.assertEqual(len(published_symbols), 18164)
        self.assertEqual(published_symbols, source_symbols)


if __name__ == "__main__":
    unittest.main()
