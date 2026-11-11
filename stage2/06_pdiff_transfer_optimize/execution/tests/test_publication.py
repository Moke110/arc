"""Non-retrieval exposure rendering survives metric consolidation."""
from pathlib import Path
import csv
import tempfile
import sys
import unittest
REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
STEP = REPO / 'stage2/06_pdiff_transfer_optimize'
sys.path.insert(0, str(STEP / 'execution/script'))
from exposure import render_exposure


class PublicationTest(unittest.TestCase):
    def test_exposure_figure_preserves_original_metadata_and_bytes(self):
        with (STEP / 'source_exposure.csv').open() as stream:
            exposure = list(csv.DictReader(stream))
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)
            render_exposure(exposure, destination)
            path = destination / '02-tahoe-exposure.png'
            self.assertEqual(path.read_bytes(), (STEP / 'output/02-tahoe-exposure.png').read_bytes())
            self.assertEqual({p.name for p in destination.iterdir()}, {'02-tahoe-exposure.png'})
