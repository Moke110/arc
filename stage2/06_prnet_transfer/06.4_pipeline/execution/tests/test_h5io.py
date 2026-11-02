"""Header-only h5ad metadata: ``open_obs_var`` must never materialise layers.

Regression for the Phase-B trap: ``anndata.read_h5ad(..., backed="r")`` keeps
``X`` lazy but eagerly loads every ``layers`` entry — 8.1 GB RSS on the frozen
Tahoe view (67,018 x 18,164) before any training starts, which is what ticket #3
hit and what ``docs/agents/execution-environment.md`` records.  The loader reads
the ``obs``/``var`` element groups with h5py instead; these tests attack that
boundary rather than trusting the implementation to stay lazy.
"""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

import _support
from ptlib import h5io


class OpenObsVarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.small_fixture_dir()

    def test_never_calls_the_backed_anndata_reader(self):
        """Poison ``anndata.read_h5ad``: the metadata path must not need it."""
        original = ad.read_h5ad

        def poison(*args, **kwargs):  # noqa: ANN002, ANN003
            raise AssertionError(
                "ptlib.h5io.open_obs_var must not call anndata.read_h5ad: "
                "backed='r' materialises every layers entry (8.1 GB on the frozen view)"
            )

        ad.read_h5ad = poison
        try:
            for name in ("cm_transfer.h5ad", "tahoe_transfer.h5ad"):
                meta = h5io.open_obs_var(self.dir / name)
                self.assertGreater(meta.n_obs, 0)
                self.assertGreater(meta.n_vars, 0)
                self.assertEqual(len(meta.obs), meta.n_obs)
                self.assertEqual(len(meta.var), meta.n_vars)
        finally:
            ad.read_h5ad = original

    def test_decodes_only_obs_and_var_elements(self):
        """X / layers are never passed to the element decoder (only obs and var)."""
        decoded: list[str | None] = []
        original = h5io._read_elem

        def spy(payload):
            decoded.append(getattr(payload, "name", None))
            return original(payload)

        h5io._read_elem = spy
        try:
            meta = h5io.open_obs_var(self.dir / "tahoe_transfer.h5ad")
        finally:
            h5io._read_elem = original
        self.assertEqual(decoded, ["/obs", "/var"])
        # ... and the header shape still comes from the X group's attrs, not values.
        adata = ad.read_h5ad(self.dir / "tahoe_transfer.h5ad")
        self.assertEqual((meta.n_obs, meta.n_vars), (adata.n_obs, adata.n_vars))

    def test_values_match_the_anndata_reader(self):
        for name in ("cm_transfer.h5ad", "tahoe_transfer.h5ad"):
            path = self.dir / name
            meta = h5io.open_obs_var(path)
            adata = ad.read_h5ad(path)  # small fixture: safe to materialise in tests
            pd.testing.assert_frame_equal(meta.obs, adata.obs, check_like=True)
            pd.testing.assert_frame_equal(meta.var, adata.var, check_like=True)
            self.assertEqual(set(meta.uns_keys), set(adata.uns.keys()))

    def test_missing_artifact_fails_closed(self):
        with self.assertRaises(h5io.H5Error):
            h5io.open_obs_var(self.dir / "does_not_exist.h5ad")

    def test_shape_disagreement_fails_closed(self):
        """A corrupted X header must fail loud, not be silently ignored."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        target = Path(tmp.name) / "tahoe_transfer.h5ad"
        shutil.copyfile(self.dir / "tahoe_transfer.h5ad", target)
        with h5py.File(target, "r+") as handle:
            shape = list(handle["X"].attrs["shape"])
            handle["X"].attrs["shape"] = [int(shape[0]) + 1, int(shape[1])]
        with self.assertRaises(h5io.H5Error):
            h5io.open_obs_var(target)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
