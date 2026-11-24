"""Stage the flat release from the accepted full payloads without rebuilding them."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
from hashlib import sha256 as hash_bytes
from pathlib import Path

# Reuse the builder's imports and release layout.
from build_views import (PAYLOAD_SCHEMA, RECIPE, SCHEMA, SOURCES,
                         open_view, release, sha256)

# Keep the original paths when regenerating the sealed #302 receipt.
PREVIOUS = Path("stage2/03_view/03.1_hvg-deg_views/output/raw-first-v1")
CURRENT = Path("stage2/03_view/03.1_hvg-deg-5k/output/release")
ARCHIVED = Path("stage2/03_view/03.1_hvg-deg-5k/historical/raw-first-v1")


def migration_manifest(previous):
    if previous.get("schema") != PAYLOAD_SCHEMA or previous.get("sources") != list(SOURCES):
        raise ValueError("expected the two-source raw-first predecessor release")
    if (previous.get("axis_size") != 5000
            or previous.get("default_family") != "full"
            or any(previous.get(field) != RECIPE[field] for field in ("normalization", "control"))):
        raise ValueError("predecessor does not use the accepted full 5k recipe")
    mapping = {
        f"{source}.{suffix}": f"{prefix}{source}.{suffix}"
        for source in SOURCES
        for suffix, prefix in (("h5ad", "full/"), ("members.csv", "full/"),
                               ("exclusions.csv", ""))
    }
    manifest = {field: previous[field] for field in
                ("axis_size", "sources", "corpus_sha256", "axis_sha256")}
    manifest.update(RECIPE, schema=SCHEMA,
                    outputs={new: previous["outputs"][old] for new, old in mapping.items()})
    return manifest, mapping


def _write_receipt(path, receipt):
    """Link a complete receipt into place without replacing an existing file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent) as staged:
        staged.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        staged.flush()
        os.link(staged.name, path)


def stage_migration(predecessor, candidate, receipt_path):
    """Create a verified candidate and a fresh receipt outside both release trees."""
    predecessor = Path(predecessor).resolve()
    candidate = Path(candidate).resolve()
    requested_receipt = Path(receipt_path)
    receipt_path = requested_receipt.resolve()
    if predecessor.is_relative_to(candidate) or candidate.is_relative_to(predecessor):
        raise ValueError("predecessor and candidate must be separate directories")
    if receipt_path.is_relative_to(predecessor) or receipt_path.is_relative_to(candidate):
        raise ValueError("migration receipt must be outside the predecessor and candidate")
    if receipt_path.exists() or requested_receipt.is_symlink():
        raise FileExistsError(f"migration receipt already exists: {requested_receipt}")

    old_manifest = (predecessor / "manifest.json").read_bytes()
    manifest, mapping = migration_manifest(json.loads(old_manifest))
    candidate.mkdir()
    try:
        for new, old in mapping.items():
            shutil.copyfile(predecessor / old, candidate / new)
        new_manifest = release.write_json(candidate / "manifest.json", manifest)
        release.check_manifest(candidate, manifest)
        for source in SOURCES:
            with open_view(candidate, source=source, corpus_sha256=manifest["corpus_sha256"],
                           axis_sha256=manifest["axis_sha256"]) as view:
                if view.shape[0] == 0:
                    raise ValueError(f"{source} has no accepted expression profiles")
                for layer in ("X", "control", "delta"):
                    view.read(0, min(view.shape[0], 64), layer)
        receipt = {
            "previous_manifest": {"path": (PREVIOUS / "manifest.json").as_posix(),
                                  "archive": (ARCHIVED / "manifest.json").as_posix(),
                                  "sha256": hash_bytes(old_manifest).hexdigest()},
            "release_manifest": {"path": (CURRENT / "manifest.json").as_posix(),
                                 "sha256": sha256(new_manifest)},
            "assets": [{"previous": (PREVIOUS / old).as_posix(),
                        "release": (CURRENT / new).as_posix(),
                        "sha256": manifest["outputs"][new]} for new, old in sorted(mapping.items())],
        }
        _write_receipt(receipt_path, receipt)
    except BaseException:
        shutil.rmtree(candidate)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("predecessor", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    stage_migration(args.predecessor, args.candidate, args.receipt)


if __name__ == "__main__":
    main()
