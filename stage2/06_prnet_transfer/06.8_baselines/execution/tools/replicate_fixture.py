#!/usr/bin/env python
"""Test aid: row-replicated copy of the synthetic frozen-input fixture.

Why: the synthetic fixture is 204 observations (204 × 18,164 cells), which is
too small to exercise the row-level allocations ticket #8 removed — the fixture
dry-run peak is dominated by imports.  The real 06.3 CM main set is ~3,700
treated QC-pass wells, so this tool repeats every fixture row ``N`` times to put
the baselines on a realistic **row count** (same gene axis, same conditions, same
fold assignment; within-condition replicates simply multiply).  It is a memory
gate input, never evidence.

``tools/fixture.py`` itself cannot be row-scaled through ``--n-main-drugs``:
above ~150 drugs per plate its ``next_well`` name generator wraps its 26 × 12
well grid and produces duplicate well ids (an anndata write error).  Replication
sidesteps that without touching ticket #4's fixture.

The h5ad is rewritten; every other frozen artifact is copied byte-for-byte and
the fixture ``manifest.json`` is re-hashed so ``ptlib.data.open_frozen_input``
verifies it.  Output stays under ``/tmp``; nothing here is tracked.

Usage::

    PY=/home/chang/miniconda3/envs/arc/bin/python
    $PY execution/tools/replicate_fixture.py \\
        --src /tmp/t08-baselines-fixture --out /tmp/t08-baselines-stress \\
        --replicates 18
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _replicated_obs(obs: pd.DataFrame, replicates: int) -> pd.DataFrame:
    """Repeat *obs* ``replicates`` times with globally unique well ids."""
    blocks = []
    for replicate in range(replicates):
        block = obs.copy()
        block["well"] = [f"{well}::rep{replicate}" for well in block["well"].astype(str)]
        block.index = pd.Index(block["well"].to_numpy(), name="well")
        blocks.append(block)
    combined = pd.concat(blocks, ignore_index=False)
    if not combined.index.is_unique:  # pragma: no cover - defensive
        raise SystemExit("replicated well ids are not unique")
    return combined


def replicate_fixture(src: Path, out: Path, replicates: int) -> dict[str, object]:
    src, out = Path(src), Path(out)
    if replicates < 1:
        raise SystemExit("--replicates must be >= 1")
    src_h5ad = src / "cm_transfer.h5ad"
    if not src_h5ad.is_file():
        raise SystemExit(f"not a frozen-input fixture: {src_h5ad} is missing")
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    adata = ad.read_h5ad(src_h5ad)
    counts = np.asarray(adata.X)
    log1p = np.asarray(adata.layers["log1p_cp10k"])
    obs = _replicated_obs(adata.obs, replicates)
    tiled_counts = np.tile(counts, (replicates, 1))
    tiled_log1p = np.tile(log1p, (replicates, 1))
    del counts, log1p

    rebuilt = ad.AnnData(X=tiled_counts, obs=obs, var=adata.var.copy())
    rebuilt.layers["log1p_cp10k"] = tiled_log1p
    for key, value in adata.uns.items():  # pragma: no cover - fixture has no uns
        rebuilt.uns[key] = value
    rebuilt.write_h5ad(out / "cm_transfer.h5ad")
    n_obs = int(rebuilt.n_obs)
    n_vars = int(rebuilt.n_vars)
    del rebuilt, tiled_counts, tiled_log1p, adata

    for path in sorted(src.iterdir()):
        if path.is_file() and path.name not in ("cm_transfer.h5ad", "manifest.json"):
            shutil.copy2(path, out / path.name)

    manifest_file = src / "manifest.json"
    payload = json.loads(manifest_file.read_text(encoding="utf-8"))
    artifacts = payload.get("artifacts", {})
    for name in sorted(artifacts):
        digest = sha256_file(out / name)
        record = artifacts[name]
        if isinstance(record, dict):
            record["sha256"] = digest
        else:
            artifacts[name] = digest
    (out / "manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        "src": str(src),
        "out": str(out),
        "replicates": int(replicates),
        "n_obs": n_obs,
        "n_vars": n_vars,
        "h5ad_bytes": int((out / "cm_transfer.h5ad").stat().st_size),
    }



def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--src", type=Path, required=True, help="source fixture directory")
    parser.add_argument("--out", type=Path, required=True, help="output fixture directory")
    parser.add_argument(
        "--replicates", type=int, default=18, help="copies of every observation"
    )
    ns = parser.parse_args(argv)
    report = replicate_fixture(ns.src, ns.out, ns.replicates)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
