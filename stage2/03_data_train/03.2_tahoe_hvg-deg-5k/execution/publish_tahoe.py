"""Expose the sealed well-level Tahoe HVG-DEG 5k data on a source-specific route."""
from __future__ import annotations

import argparse
from errno import EXDEV
import json
from pathlib import Path
import shutil
import sys

import numpy as np

STEP = Path(__file__).resolve().parents[1]
TRAINING = STEP.parent
STAGE = TRAINING.parent
ARCHIVE = TRAINING / "03.1_hvg-deg-5k/output"
DEFAULT = STEP / "output/release"
STAGING = STAGE.parent / "tmp/stage2-tahoe-training-staging"
SOURCE = "tahoe_bulk"
sys.path[:0] = [str(ARCHIVE / "api"), str(STAGE / "util/output/api")]

from expression_view import RECIPE, SCHEMA, check_inventory, open_view, output_names
from stage2_util.derivation import release
from stage2_util.publication import publication_lock, publish_candidate


def _manifest(sealed):
    path = sealed / "manifest.json"
    manifest = release.read_json(path)
    if (manifest.get("schema") != SCHEMA or set(check_inventory(manifest)) != {"tenv", SOURCE}
            or any(manifest.get(key) != value for key, value in RECIPE.items())):
        raise ValueError("expected the sealed combined full-view release")
    return {
        **manifest, "sources": [SOURCE], "sealed_manifest_sha256": release.sha256(path),
        "outputs": {name: manifest["outputs"][name] for name in output_names((SOURCE,))},
    }


def verify_tahoe(sealed, destination):
    """Verify the unchanged Tahoe well payload and both ledgers against the sealed manifest."""
    destination = Path(destination)
    expected = _manifest(Path(sealed))
    path = destination / "manifest.json"
    if (release.read_json(path) != expected
            or {asset.name for asset in destination.iterdir()} != set(expected["outputs"]) | {"manifest.json"}):
        raise ValueError("Tahoe release differs from the sealed well-level inventory")
    name = f"{SOURCE}.members.csv"
    if release.sha256(destination / name) != expected["outputs"][name]:
        raise ValueError("Tahoe membership ledger digest differs from the sealed release")
    with open_view(destination, source=SOURCE, corpus_sha256=expected["corpus_sha256"],
                   axis_sha256=expected["axis_sha256"]) as view:
        stop = min(64, view.shape[0])
        if not np.array_equal(view.read(0, stop, "delta"), view.read(0, stop) - view.read(0, stop, "control")):
            raise ValueError("Tahoe well-level Delta is inconsistent")
        return {"result": "PASS", "manifest_sha256": view.manifest_sha256,
                "shape": list(view.shape), "outputs": expected["outputs"]}


def publish_tahoe(sealed, destination, *, staging=STAGING):
    """Atomically publish a digest-equal alias of the archived Tahoe well-level assets."""
    sealed = Path(sealed).resolve()
    expected = _manifest(sealed)

    def build(candidate):
        for name in expected["outputs"]:
            source, target = sealed / name, candidate / name
            try:
                target.hardlink_to(source)
            except OSError as error:
                if error.errno != EXDEV:
                    raise
                shutil.copyfile(source, target)
        release.write_json(candidate / "manifest.json", expected)

    with publication_lock(destination) as output:
        stage = Path(staging).resolve()
        if (sealed.is_relative_to(output) or output.is_relative_to(sealed)
                or stage.is_relative_to(output) or stage.is_relative_to(sealed)):
            raise ValueError("Tahoe archive, destination and staging must be separate")
        if output.exists() and not (output / "manifest.json").is_file():
            raise ValueError("refusing to replace unrecognized Tahoe output")
        return publish_candidate(release, output, build, lambda root: verify_tahoe(sealed, root),
                                 staging=stage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("publish", "verify"))
    parser.add_argument("--sealed", type=Path, default=ARCHIVE / "release")
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--staging", type=Path, default=STAGING)
    args = parser.parse_args()
    if args.output.resolve() == DEFAULT.resolve() and release.sha256(args.sealed / "manifest.json") != (
        release.sha256(ARCHIVE / "release/manifest.json")
    ):
        parser.error("stable Tahoe publication requires the accepted full-view manifest")
    if args.command == "publish":
        output = publish_tahoe(args.sealed, args.output, staging=args.staging)
        result = {"result": "PASS", "manifest_sha256": release.sha256(output / "manifest.json")}
    else:
        result = verify_tahoe(args.sealed, args.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
