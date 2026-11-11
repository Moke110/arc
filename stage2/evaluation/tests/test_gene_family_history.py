"""Label preparation keeps its original gene family after Round 1 replacement."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from stage2.evaluation.build_labels import significance_gene_family


class GeneFamilyHistoryTest(unittest.TestCase):
    def test_active_then_archived_family_has_identical_order(self):
        with TemporaryDirectory() as temp:
            step=Path(temp)
            active=step/'output/round1/labels/test/genes.csv'
            historical=step/'execution/history/round1-exploratory/round1/labels/test/genes.csv'
            active.parent.mkdir(parents=True)
            active.write_text('symbol,hgnc_id\nA,H1\nHDGFL3,H2\n')
            original=significance_gene_family(step)
            historical.parent.mkdir(parents=True)
            active.rename(historical)
            self.assertTrue(original.equals(significance_gene_family(step)))
            historical.unlink()
            with self.assertRaises(FileNotFoundError):
                significance_gene_family(step)
