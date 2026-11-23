"""Stage the flat release from the accepted full payloads without rebuilding them."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

# Reuse the builder's imports and release layout.
from build_views import (PAYLOAD_SCHEMA, RECIPE, SCHEMA, SOURCES,
                         open_view, release, sha256)

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


def stage_migration(predecessor, candidate):
    """Create and verify a new directory, leaving the predecessor available."""
    predecessor, candidate = Path(predecessor), Path(candidate)
    old_manifest = predecessor / "manifest.json"
    manifest, mapping = migration_manifest(release.read_json(old_manifest))
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
        return {
            "previous_manifest": {"path": (PREVIOUS / "manifest.json").as_posix(),
                                  "archive": (ARCHIVED / "manifest.json").as_posix(),
                                  "sha256": sha256(old_manifest)},
            "release_manifest": {"path": (CURRENT / "manifest.json").as_posix(),
                                 "sha256": sha256(new_manifest)},
            "assets": [{"previous": (PREVIOUS / old).as_posix(),
                        "release": (CURRENT / new).as_posix(),
                        "sha256": manifest["outputs"][new]} for new, old in sorted(mapping.items())],
        }
    except BaseException:
        shutil.rmtree(candidate)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("predecessor", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    release.write_json(args.receipt, stage_migration(args.predecessor, args.candidate))


if __name__ == "__main__":
    main()
