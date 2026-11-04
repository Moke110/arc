"""G2/G3 contract for the seven crosswalk tables and the coverage README (#16).

Run from the repository root:

    python -m unittest discover -s stage2/02_dataset_corpus/execution/tests -v

The tables and README are small tracked artifacts, so the inventory tests run
even when the local-only h5ad corpora are absent; the coverage tests need the
five corpora (and skip otherwise).  Every assertion reads CSV/JSON/Markdown or
an h5ad ``var``/``obs`` header — never a corpus matrix.
"""

from __future__ import annotations

import hashlib
import json
import sys
import unittest
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from cmdata import paths  # noqa: E402

TABLES = paths.TABLES_DIR
README = paths.CORPUS_V2_README
RECORD = paths.COVERAGE_RECORD
TAHOE = paths.TAHOE_H5AD
CM_CORPORA = {corpus_id: paths.corpus_h5ad(corpus_id)
              for corpus_id in paths.CORPUS_IDS}
ALL_CORPORA = [*paths.CORPUS_IDS, "tahoe"]
FROZEN = paths.REPO_ROOT / "stage2" / "_frozen" / "cm-data-v1"
CROSSWALK_FOUR_SOURCE_ROWS = 143647          # frozen v1 crosswalk row count
CROSSWALK_TAHOE_ROWS = 125420                # 2 native id namespaces x 62,710
CORPUS_COVERAGE_COLUMNS = [
    "corpus", "distinct cell lines", "chemicals (perturbed wells)",
    "distinct doses (µM)", "time window", "(cell line, chemical, dose)",
    "control wells", "blank wells",
]
T24_COLUMNS = [
    "corpus", "distinct cell lines", "chemicals", "distinct doses",
    "(cell line, chemical, dose)", "basis",
]


def _missing(paths_: list[Path]) -> list[Path]:
    return [path for path in paths_ if not path.exists()]


def _tables_skip() -> str:
    missing = _missing([TABLES / name for name in paths.TABLE_FILES]
                       + [README, RECORD])
    return "" if not missing else f"tables/coverage not built: {missing}"


def _corpora_skip() -> str:
    missing = _missing([*CM_CORPORA.values(), TAHOE])
    return "" if not missing else f"corpus h5ads not present locally: {missing}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_crosswalk() -> pd.DataFrame:
    return pd.read_csv(TABLES / "gene_crosswalk.csv", dtype=str,
                       keep_default_na=False)


def _read_readme_table(header: list[str]) -> list[dict[str, str]]:
    """Parse the first Markdown table whose header matches *header*."""
    lines = README.read_text().splitlines()
    header_line = "| " + " | ".join(header) + " |"
    rows: list[dict[str, str]] = []
    found = False
    for line in lines:
        if not found:
            found = line.strip() == header_line
            continue
        if not line.startswith("|"):
            break
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != len(header):
            continue
        if all(set(cell) <= {"-", ":"} for cell in cells):
            continue
        rows.append(dict(zip(header, cells)))
    if not found:
        raise AssertionError(f"README table not found: {header_line}")
    return rows


def _number(cell: str) -> int:
    return int(cell.replace(",", ""))


@unittest.skipUnless(not _tables_skip(), _tables_skip())
class TablesInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads(RECORD.read_text())

    def test_exactly_seven_tables(self):
        present = sorted(path.name for path in TABLES.iterdir()
                         if path.is_file())
        self.assertEqual(present, sorted(paths.TABLE_FILES))

    def test_row_counts(self):
        genes = pd.read_csv(TABLES / "genes.csv", dtype=str)
        self.assertEqual(len(genes), 18164)
        self.assertEqual(genes["symbol"].nunique(), 18164)
        compounds = pd.read_csv(TABLES / "compounds.csv",
                                keep_default_na=False)
        self.assertEqual(len(compounds), 532)
        self.assertEqual(compounds["drug_key"].nunique(), 525)
        excluded = pd.read_csv(TABLES / "excluded_compounds.csv",
                               keep_default_na=False)
        non_blank = excluded[excluded["entity_type"] != "technical_blank"]
        blanks = excluded[excluded["entity_type"] == "technical_blank"]
        self.assertEqual(len(non_blank), 8)
        self.assertEqual(int(non_blank["n_wells"].astype(int).sum()), 143)
        self.assertEqual(len(blanks), 1)
        self.assertEqual(int(blanks["n_wells"].astype(int).sum()), 32)
        source_sets = pd.read_csv(TABLES / "source_gene_sets.csv", dtype=str)
        self.assertEqual(len(source_sets), 104149)
        metadata_cm = pd.read_csv(TABLES / "metadata_cm.csv",
                                  keep_default_na=False)
        self.assertEqual(len(metadata_cm), 8673)
        metadata_tahoe = pd.read_csv(TABLES / "metadata_tahoe.csv",
                                     keep_default_na=False)
        self.assertEqual(len(metadata_tahoe), 67018)
        crosswalk = _read_crosswalk()
        self.assertEqual(len(crosswalk), CROSSWALK_FOUR_SOURCE_ROWS
                         + CROSSWALK_TAHOE_ROWS)
        self.assertEqual(int((crosswalk["dataset_id"] != "tahoe").sum()),
                         CROSSWALK_FOUR_SOURCE_ROWS)
        self.assertEqual(int((crosswalk["dataset_id"] == "tahoe").sum()),
                         CROSSWALK_TAHOE_ROWS)

    def test_metadata_columns_and_ids(self):
        metadata_cm = pd.read_csv(TABLES / "metadata_cm.csv",
                                  keep_default_na=False)
        expected = ["well", "dataset_id", "legacy_dataset_id", "source_file",
                    "accession", *paths.CORPUS_OBS[2:]]
        self.assertEqual(list(metadata_cm.columns), expected)
        self.assertEqual(set(metadata_cm["dataset_id"]), set(paths.CORPUS_IDS))
        legacy = {"dtoxs_primary": "dtoxs_cm", "dtoxs_ipsc": "dtoxs_cm",
                  "sctox": "screen_ctox", "tenv": "temp_env"}
        for corpus_id, legacy_id in legacy.items():
            rows = metadata_cm[metadata_cm["dataset_id"] == corpus_id]
            self.assertEqual(len(rows), paths.EXPECTED_CORPUS_COUNTS[corpus_id]
                             ["rows"])
            self.assertEqual(set(rows["legacy_dataset_id"]), {legacy_id})
            self.assertTrue((rows["source_file"].astype(str) != "").all())
            self.assertTrue((rows["accession"].astype(str) != "").all())
        self.assertTrue(metadata_cm["well"].is_unique)
        metadata_tahoe = pd.read_csv(TABLES / "metadata_tahoe.csv",
                                     keep_default_na=False)
        expected = ["profile_id", "dataset_id", "is_blank", *paths.TAHOE_OBS]
        self.assertEqual(list(metadata_tahoe.columns), expected)
        self.assertEqual(set(metadata_tahoe["dataset_id"]), {"tahoe"})
        self.assertFalse(
            (metadata_tahoe["is_blank"].astype(str) != "False").any())
        self.assertTrue(metadata_tahoe["profile_id"].is_unique)

    def test_file_and_content_hashes_match_record(self):
        self.assertEqual(_sha256(README),
                         self.record["readme"]["sha256"])
        for name in paths.TABLE_FILES:
            info = self.record["tables"][name]
            self.assertEqual(_sha256(TABLES / name), info["sha256"])
            self.assertEqual((TABLES / name).stat().st_size, info["bytes"])
            self.assertEqual(info["rows"],
                             len(pd.read_csv(TABLES / name,
                                             keep_default_na=False)))

    def test_frozen_identity_tables_are_content_identical(self):
        checks = self.record["frozen_identity_check"]
        if checks.get("genes.csv") is None:
            self.skipTest("frozen v1 CM evidence is not present")
        self.assertTrue(checks["genes.csv"])
        self.assertTrue(checks["compounds.csv"])
        self.assertTrue(checks["source_gene_sets.csv"])
        self.assertTrue(checks["gene_crosswalk.csv (four-source part)"])
        frozen = FROZEN / "corpus" / "v1"
        self.assertEqual((TABLES / "genes.csv").read_bytes(),
                         (frozen / "genes.csv").read_bytes())
        self.assertEqual((TABLES / "compounds.csv").read_bytes(),
                         (frozen / "compounds.csv").read_bytes())
        self.assertEqual((TABLES / "source_gene_sets.csv").read_bytes(),
                         (frozen / "source_gene_sets.csv").read_bytes())
        crosswalk = _read_crosswalk()
        four_source = crosswalk[crosswalk["dataset_id"] != "tahoe"]
        self.assertEqual(
            four_source.to_csv(index=False, lineterminator="\n").encode(),
            (FROZEN / "datasets-1" / "gene_crosswalk.csv").read_bytes())

    def test_readme_coverage_matches_record(self):
        coverage = self.record["corpora"]
        rows = {row["corpus"].strip("`"): row
                for row in _read_readme_table(CORPUS_COVERAGE_COLUMNS)}
        self.assertEqual(set(rows), set(ALL_CORPORA))
        for corpus_id in ALL_CORPORA:
            entry = coverage[corpus_id]
            row = rows[corpus_id]
            self.assertEqual(_number(row["distinct cell lines"]),
                             entry["distinct_cell_lines"])
            self.assertEqual(_number(row["chemicals (perturbed wells)"]),
                             entry["chemicals"])
            self.assertEqual(_number(row["distinct doses (µM)"]),
                             entry["distinct_doses"])
            self.assertEqual(_number(row["(cell line, chemical, dose)"]),
                             entry["combinations"])
            self.assertEqual(row["time window"], entry["time_window"])
            self.assertEqual(_number(row["control wells"]), entry["control"])
            self.assertEqual(_number(row["blank wells"]), entry["blank"])

    def test_readme_24h_subrow(self):
        coverage = self.record["corpora"]
        rows = {row["corpus"].strip("`"): row
                for row in _read_readme_table(T24_COLUMNS)}
        self.assertEqual(set(rows), set(ALL_CORPORA))
        for corpus_id in ("dtoxs_primary", "dtoxs_ipsc"):
            self.assertIsNone(coverage[corpus_id]["t24"])
            for column in T24_COLUMNS[1:-1]:
                self.assertEqual(rows[corpus_id][column], "n.a.")
            self.assertIn("n.a.", rows[corpus_id]["basis"])
        for corpus_id in ("sctox", "tenv", "tahoe"):
            entry = coverage[corpus_id]["t24"]
            self.assertIsNotNone(entry)
            row = rows[corpus_id]
            self.assertEqual(_number(row["distinct cell lines"]),
                             entry["distinct_cell_lines"])
            self.assertEqual(_number(row["chemicals"]), entry["chemicals"])
            self.assertEqual(_number(row["distinct doses"]),
                             entry["distinct_doses"])
            self.assertEqual(_number(row["(cell line, chemical, dose)"]),
                             entry["combinations"])
            self.assertEqual(row["basis"], entry["basis"])

    def test_readme_rename_table(self):
        text = README.read_text()
        for token in ("`screen_ctox`", "`sctox`", "`temp_env`", "`tenv`",
                      "`dtoxs_cm`", "`dtoxs_primary`", "`dtoxs_ipsc`",
                      "`04_tahoe_data`", "`02_dataset_corpus`"):
            self.assertIn(token, text)
        self.assertIn("Historical records keep the **old** names", text)

    def test_coverage_record_shape(self):
        self.assertEqual(self.record["ticket"], "#16")
        self.assertEqual(set(self.record["corpora"]), set(ALL_CORPORA))
        for corpus_id in ALL_CORPORA:
            entry = self.record["corpora"][corpus_id]
            for key in ("rows", "perturbed", "control", "blank",
                        "distinct_cell_lines", "chemicals", "distinct_doses",
                        "time_window", "combinations"):
                self.assertIn(key, entry)
        for corpus_id in ("dtoxs_primary", "dtoxs_ipsc"):
            self.assertIsNone(self.record["corpora"][corpus_id]["t24"])
        # Ticket #17 replaced the "G4-G6 not run" placeholder with the
        # executable gates' results and the verify record pointer.
        self.assertNotIn("G4-G6", self.record["gate_results"])
        for gate in ("G4", "G5", "G6"):
            self.assertTrue(self.record["gate_results"][gate]["passed"])
        self.assertEqual(self.record["gate_results"]["verified_by"],
                         "stage2/02_dataset_corpus/execution/build_record/"
                         "verify.json")
        self.assertEqual(self.record["verification"]["result"], "PASS")


@unittest.skipUnless(not _tables_skip() and not _corpora_skip(),
                     _tables_skip() or _corpora_skip())
class CoverageContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.crosswalk = _read_crosswalk()

    def test_crosswalk_covers_every_cm_var_id(self):
        genes = set(pd.read_csv(TABLES / "genes.csv", dtype=str)["symbol"])
        self.assertEqual(len(genes), 18164)
        for corpus_id, path in CM_CORPORA.items():
            source_id = paths.CORPUS_TO_SOURCE[corpus_id]
            with h5py.File(path, "r") as handle:
                var = read_elem(handle["var"])
            rows = self.crosswalk[self.crosswalk["dataset_id"] == source_id]
            symbols = set(var.index.astype(str))
            self.assertTrue(symbols <= set(rows["symbol"]),
                            f"{corpus_id}: var symbols missing from crosswalk")
            self.assertTrue(set(var["hgnc_id"].astype(str))
                            <= set(rows["hgnc_id"]),
                            f"{corpus_id}: var hgnc_ids missing from crosswalk")
            self.assertTrue(genes <= symbols,
                            f"{corpus_id}: 18,164 CM symbols not a subset of var")

    def test_crosswalk_covers_tahoe_native_ids(self):
        with h5py.File(TAHOE, "r") as handle:
            var = read_elem(handle["var"])
        rows = self.crosswalk[self.crosswalk["dataset_id"] == "tahoe"]
        symbols = set(var.index.astype(str))
        self.assertTrue(symbols <= set(rows["symbol"]))
        self.assertTrue(set(var["ensembl_id"].astype(str))
                        <= set(rows["source_gene_id"]))
        self.assertTrue(set(var["token_id"].astype(str))
                        <= set(rows["source_gene_id"]))
        self.assertEqual(len(var), 62710)

    def test_compounds_cover_every_perturbed_drug_key(self):
        compounds = set(pd.read_csv(
            TABLES / "compounds.csv", keep_default_na=False)["drug_key"])
        metadata = pd.read_csv(TABLES / "metadata_cm.csv",
                               keep_default_na=False)
        perturbed = metadata[(metadata["is_control"].astype(str) == "False")
                             & (metadata["is_blank"].astype(str) == "False")]
        self.assertTrue(set(perturbed["drug_key"]) <= compounds,
                        "compounds.csv does not cover every perturbed drug_key")
        # ... and Tahoe's chemical count is its own native vocabulary.
        tahoe = pd.read_csv(TABLES / "metadata_tahoe.csv",
                            keep_default_na=False)
        control = tahoe["is_control"].astype(str) == "True"
        self.assertEqual(tahoe.loc[~control, "drug"].nunique(), 379)

    def test_metadata_rows_are_the_corpora(self):
        metadata = pd.read_csv(TABLES / "metadata_cm.csv",
                               keep_default_na=False)
        for corpus_id, path in CM_CORPORA.items():
            with h5py.File(path, "r") as handle:
                obs = read_elem(handle["obs"])
            rows = metadata[metadata["dataset_id"] == corpus_id]
            self.assertEqual(len(rows), len(obs))
            self.assertEqual(set(rows["well"]), set(obs.index.astype(str)))
        with h5py.File(TAHOE, "r") as handle:
            obs = read_elem(handle["obs"])
        tahoe = pd.read_csv(TABLES / "metadata_tahoe.csv",
                            keep_default_na=False)
        self.assertEqual(len(tahoe), len(obs))
        self.assertEqual(set(tahoe["profile_id"]), set(obs.index.astype(str)))


if __name__ == "__main__":
    unittest.main()
