"""Acceptance checks at the Denoise-DEG evidence and release boundaries."""

from pathlib import Path
import sys

import pandas as pd

STEP = Path(__file__).resolve().parents[2]
REPO = STEP.parents[2]
sys.path.insert(0, str(STEP / "execution"))

from diagnose_scores import disease_evidence


def test_votes_use_complete_eligible_published_deg1000_lists():
    disease, _ = disease_evidence(REPO)
    assert disease.groupby("disease", sort=False).size().to_dict() == {
        "nHF": 1000,
        "ICM": 1000,
        "genetic_cardiomyopathy": 1000,
        "TOX": 1000,
    }
    assert disease.loc[disease.disease.eq("nHF")].iloc[0].gene_index == 16837
    assert not disease.duplicated(["disease", "gene_index"]).any()


def test_release_reports_full_list_coverage_and_frozen_nested_axes():
    reference, _ = disease_evidence(REPO)
    report = pd.read_csv(STEP / "figure/data/disease-DEG1000-coverage.csv")
    assert set(report.axis_size) == {1000, 2000, 3000, 4000, 5000}
    previous = set()
    for size, rows in report.groupby("axis_size"):
        axis = pd.read_csv(STEP / f"output/release/gene_axis_{size}.csv")
        assert len(axis) == axis.gene_index.nunique() == size
        assert axis.gene_index.is_monotonic_increasing
        assert rows.disease.nunique() == len(rows) == 4
        selected = set(axis.gene_index)
        assert previous < selected
        for row in rows.itertuples(index=False):
            actual = set(reference.loc[reference.disease.eq(row.disease), "gene_index"])
            assert row.total == len(actual)
            assert row.hits == len(actual & selected)
            assert row.coverage == row.hits / row.total
        previous = selected
