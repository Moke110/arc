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
REPO = next(parent for parent in STEP.parents if (parent / ".git").exists())
STAGE = REPO / "cpst"
DEFAULT = STEP / "output/release"
STAGING = REPO / "tmp/stage2-tahoe-training-staging"
sys.path[:0] = [str(STEP / "output/api"), str(STAGE / "util/output/api")]

from tahoe_expression import (ACCEPTED_TAHOE_SHA256, RECIPE, SCHEMA, SOURCE,
                              check_inventory, open_tahoe, output_names)
from cpst_util.derivation import release
from cpst_util.publication import publication_lock, publish_candidate


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


def verify_release(destination, *, expected_sha256=ACCEPTED_TAHOE_SHA256):
    """Verify the accepted Tahoe release without reading the combined view."""
    destination = Path(destination)
    with open_tahoe(destination, expected_sha256=expected_sha256) as view:
        outputs = view.manifest["outputs"]
        if {asset.name for asset in destination.iterdir()} != set(outputs) | {"manifest.json"}:
            raise ValueError("Tahoe release differs from the sealed well-level inventory")
        name = f"{SOURCE}.members.csv"
        if release.sha256(destination / name) != outputs[name]:
            raise ValueError("Tahoe membership ledger digest differs from the sealed release")
        stop = min(64, view.shape[0])
        if not np.array_equal(view.read(0, stop, "delta"), view.read(0, stop) - view.read(0, stop, "control")):
            raise ValueError("Tahoe well-level Delta is inconsistent")
        return {"result": "PASS", "manifest_sha256": view.manifest_sha256,
                "shape": list(view.shape), "outputs": outputs}


def verify_tahoe(sealed, destination):
    """Compare a reconstructed release against its explicit sealed source."""
    destination = Path(destination)
    expected = _manifest(Path(sealed))
    if release.read_json(destination / "manifest.json") != expected:
        raise ValueError("Tahoe release differs from the sealed well-level inventory")
    return verify_release(destination, expected_sha256=None)


def publish_tahoe(sealed, destination, *, staging=STAGING):
    """Atomically publish a digest-equal alias of the archived Tahoe well-level assets."""
    sealed = Path(sealed).resolve()
    expected = _manifest(sealed)

    def verify(candidate):
        if Path(destination).resolve() == DEFAULT.resolve():
            return verify_release(candidate)
        return verify_tahoe(sealed, candidate)

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
        return publish_candidate(release, output, build, verify,
                                 staging=stage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("publish", "verify"))
    parser.add_argument("--sealed", type=Path, help="combined full view for historical reconstruction")
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--staging", type=Path, default=STAGING)
    args = parser.parse_args()
    if args.command == "publish":
        if args.sealed is None:
            parser.error("publish requires --sealed")
        output = publish_tahoe(args.sealed, args.output, staging=args.staging)
        result = {"result": "PASS", "manifest_sha256": release.sha256(output / "manifest.json")}
    else:
        result = (verify_release(args.output) if args.sealed is None
                  else verify_tahoe(args.sealed, args.output))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
