"""Step-level acceptance invariants across the six released disease-deg catalogs.

This is #192's cross-category gate. The per-category tests (#187, #188, #189) each
check their own catalog; this file checks the property that only becomes visible
when all six are read together: one canonical ``stage2-disease-deg-catalog-v1``
shape, a self-consistent ``catalog_sha256``, a declared case/control basis, and
digest-pinned ``derived_inputs`` where a catalog declares them.

**Documented cross-ticket shape variance (accepted, not silently tolerated).**
The six catalogs were produced by two tickets and are *not* byte-for-byte the same
shape. The test pins the variance explicitly, so a new third shape fails the gate:

* ``af``/``icm``/``nhf``/``genetic_cardiomyopathy`` (#187/#188) declare
  ``axis: {path: output/gene-identity/genes.csv, ...}`` (the published copy with
  ``gene_index``), a ``rank_contract`` block, a ``rejected`` list, ``inputs``
  values as **bare 64-hex sha256 strings**, and a ``derived_inputs`` list of
  ``{relative, sha256}`` records for gitignored item-local artefacts.
* ``cardiotox_TOP2``/``cardiotox_other`` (#189) declare
  ``axis: {path: output/gene-identity/tables/genes.csv, ...}`` (the source table),
  no ``rank_contract``/``rejected`` keys, ``inputs`` values as
  ``{bytes, sha256}`` dicts, and ``derived_inputs: null`` (their inputs are the
  tracked item DEG JSON tables, not gitignored item-local CSVs).

Both axis files resolve the same frozen 18,164-gene axis in the same order; the
test asserts that equivalence instead of pretending the paths are identical. The
``n_case``/``n_control`` rule also admits an explicit ``*_basis`` string: #189's
pooled class units record why an arm count is what it is, and a zero/absent arm
count without a basis is a failure, never a silently tolerated null.
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
    "cardiotox_TOP2",
    "cardiotox_other",
)
TICKET_187_188 = {"af", "icm", "nhf", "genetic_cardiomyopathy"}
TICKET_189 = {"cardiotox_TOP2", "cardiotox_other"}
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

    def test_all_six_categories_have_catalogs_and_nothing_else(self):
        self.assertEqual(set(self.catalogs), set(EXPECTED_CATEGORIES))
        self.assertEqual(len(self.catalogs), len(EXPECTED_CATEGORIES))
        self.assertTrue((DISEASE_DEG / "cardiotox" / "crosswalk.csv").is_file(),
                        "the cardiotox support dir must stay catalogs-free")

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
        # ADR 0068: an all-tied dataset is published with informative: false and
        # a recorded reason; the field is a boolean, never inferred downstream.
        for category, catalog in self.catalogs.items():
            for record in catalog["datasets"]:
                with self.subTest(category=category, dataset=record["dataset_id"]):
                    self.assertIsInstance(record["informative"], bool)
                    details = record["informativeness"]
                    self.assertIsInstance(details["rule"], str)
                    self.assertGreaterEqual(int(details["n_distinct_effect"]), 1)
                    if record["informative"]:
                        self.assertIsNone(details["reason"])
                        self.assertGreater(details["n_distinct_effect"], 1)
                    else:
                        self.assertIsInstance(details["reason"], str)
                        self.assertTrue(details["reason"])
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
        """The #187/#188 bare-digest and #189 {bytes, sha256} shapes, pinned."""
        for category, catalog in self.catalogs.items():
            expected_shape = "bare" if category in TICKET_187_188 else "dict"
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
                    if category in TICKET_187_188:
                        self.assertIsInstance(derived, list)
                        self.assertGreater(len(derived), 0)
                        for entry in derived:
                            self.assertEqual(set(entry), {"relative", "sha256"})
                            self.assertFalse(Path(entry["relative"]).is_absolute())
                            self.assertTrue(HEX64.match(entry["sha256"]),
                                            f"{record['dataset_id']}: bad digest")
                    else:
                        # #189's class units pin tracked item JSON tables in
                        # `inputs`; there is no gitignored item-local artefact.
                        self.assertIsNone(derived)

    def test_axis_declarations_resolve_to_the_same_frozen_axis(self):
        """#187/#188 pin the published copy, #189 the source table; both agree."""
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
                expected = (published_digest if category in TICKET_187_188
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
