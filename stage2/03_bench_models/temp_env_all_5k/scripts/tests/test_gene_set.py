"""Offline contract for the 5,000-gene HVG selection (ticket 01).

Re-derives the selection from the persisted per-gene ranks and checks the
rank rule, so a hand-edited or non-deterministic CSV fails closed.  No h5ad is
read here.
"""

from __future__ import annotations

import hashlib
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
BENCH = HERE.parents[2]
DATA = BENCH / "data"
ROOT = BENCH.parents[2]
sys.path.insert(0, str(HERE.parents[1]))

import build_gene_set  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class GeneSetContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads((DATA / "gene_set_5k.manifest.json").read_text())
        cls.selected = pd.read_csv(DATA / "gene_set_5k.csv")
        cls.ranks = pd.read_csv(DATA / "gene_set_5k_ranks.csv")
        cls.corpus = pd.read_csv(
            ROOT / "stage2/02_cm_data/corpus/v1/genes.csv", usecols=["symbol"]
        )["symbol"].astype(str).to_numpy()

    def test_manifest_hashes(self) -> None:
        self.assertEqual(
            sha256_file(DATA / "gene_set_5k.csv"),
            self.manifest["selection"]["sha256"],
        )
        self.assertEqual(
            sha256_file(DATA / "gene_set_5k_ranks.csv"),
            self.manifest["selection"]["all_ranks"]["sha256"],
        )

    def test_selection_is_5000_unique_corpus_genes(self) -> None:
        s = self.selected
        self.assertEqual(len(s), 5000)
        self.assertTrue(s["symbol"].is_unique)
        self.assertTrue(np.all(np.diff(s["corpus_index"].to_numpy()) > 0))
        self.assertTrue(s["corpus_index"].between(0, len(self.corpus) - 1).all())
        self.assertTrue(np.array_equal(s["symbol"].to_numpy(), self.corpus[s["corpus_index"]]))

    def test_rank_sum_is_the_sum_of_three_ranks(self) -> None:
        expected = sum(self.ranks[f"{name}__rank"] for name in build_gene_set.SOURCES)
        self.assertTrue(np.array_equal(expected.to_numpy(), self.ranks["rank_sum"].to_numpy()))

    def test_non_evaluable_genes_get_the_worst_rank(self) -> None:
        for name in build_gene_set.SOURCES:
            eligible = self.ranks[f"{name}__eligible"].to_numpy(dtype=bool)
            rank = self.ranks[f"{name}__rank"].to_numpy()
            self.assertTrue((rank[~eligible] == build_gene_set.WORST_RANK).all())
            self.assertEqual(int(eligible.sum()), int(rank[eligible].max()))
            self.assertEqual(int(rank[eligible].min()), 1)

    def test_ranks_follow_descending_dispersion(self) -> None:
        for name in build_gene_set.SOURCES:
            eligible = self.ranks[self.ranks[f"{name}__eligible"].astype(bool)]
            ordered = eligible.sort_values(f"{name}__rank")
            dispersion = ordered[f"{name}__dispersions_norm"].to_numpy()
            self.assertTrue(np.all(np.diff(dispersion) <= 1e-12), name)

    def test_selected_are_the_5000_lowest_rank_sums(self) -> None:
        corpus_index = self.ranks["corpus_index"].to_numpy()
        rank_sum = self.ranks["rank_sum"].to_numpy()
        order = np.lexsort((corpus_index, rank_sum))
        expected = np.sort(corpus_index[order[:5000]])
        self.assertTrue(np.array_equal(expected, self.selected["corpus_index"].to_numpy()))

    def test_manifest_declares_no_batch_key_and_all_cells(self) -> None:
        method = self.manifest["method"]
        self.assertIsNone(method["batch_key"])
        self.assertEqual(method["flavor"], "seurat")
        self.assertEqual(method["layer"], "log1p_cp10k")
        self.assertIn("all wells", method["cells"])
        self.assertEqual(self.manifest["target_size"], 5000)


if __name__ == "__main__":
    unittest.main()
