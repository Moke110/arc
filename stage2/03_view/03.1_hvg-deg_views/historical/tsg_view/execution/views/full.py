"""The `full` view: a digest-equal rebuild of the accepted 01 release, aliased.

The rebuild runs in this step's staging directory through the owning step's own
``bases.view.build`` so the ADR 0077 recipe is never re-derived, then compares the
rebuilt manifest dictionary with the accepted manifest and re-hashes every staged byte
against the accepted digests. Nothing is published and nothing is written inside
``stage2/01_bases``; the product carries only the accepted digests and recipe.
"""
import shutil

from . import (CONTROL_RULE, RECIPE, HarnessError, relative, sha256)


def build(context):
    view = context.modules["view"]
    release = context.modules["release"]
    accepted_root = context.accepted_dir
    accepted = context.accepted_manifest()
    accepted_sha = sha256(accepted_root / "manifest.json")
    # The alias claims these bytes are authoritative, so re-hash them before rebuilding.
    try:
        release.check_manifest(accepted_root, accepted)
    except release.ReleaseError as error:
        raise HarnessError(f"the accepted release is not intact: {error}") from error

    stage = context.staging_dir / "rebuild"
    if stage.exists():
        shutil.rmtree(stage)
    rebuilt = view.build(stage, corpus_dir=context.corpus_dir, axis_dir=context.axis_dir,
                         size=context.size, name=context.name,
                         repo_root=context.provenance_root, fixture=context.fixture)
    if rebuilt != accepted:
        differing = [key for key in sorted(set(rebuilt) | set(accepted))
                     if rebuilt.get(key) != accepted.get(key)]
        raise HarnessError(
            "full-view rebuild differs from the accepted release in "
            f"{differing} (artifact: manifest.json)")
    try:
        release.check_manifest(stage, accepted)
    except release.ReleaseError as error:
        raise HarnessError(f"the rebuilt release is not digest-equal: {error}") from error
    rebuilt_sha = sha256(stage / "manifest.json")
    if rebuilt_sha != accepted_sha:
        raise HarnessError(
            "full-view rebuild manifest bytes differ from the accepted release "
            f"({rebuilt_sha} != {accepted_sha}; artifact: manifest.json)")
    context.evidence.update({
        "accepted_release": relative(accepted_root, context.root),
        "accepted_manifest_sha256": accepted_sha,
        "rebuilt_manifest_sha256": rebuilt_sha,
        "rebuilt_outputs": {name: record["sha256"]
                            for name, record in sorted(accepted["outputs"].items())},
        "digest_equal": True,
    })
    shutil.rmtree(stage, ignore_errors=True)

    with context.accepted_source() as source:
        rows = source.obs.copy()
    rows["profile_id"] = rows.index.astype(str)
    alias = {
        "release": relative(accepted_root, context.root),
        "manifest_sha256": accepted_sha,
        "payload_sha256": {name: accepted["outputs"][f"{name}.h5ad"]["sha256"]
                           for name in ("tenv", "tahoe_bulk")},
        "outputs": {name: record["sha256"]
                    for name, record in sorted(accepted["outputs"].items())},
        "recipe": accepted["recipe"],
        "normalization": accepted["normalization"],
        "control": accepted["control"],
    }
    return context.write_product(rows=rows, blocks=None, mode="alias",
                                 recipe=RECIPE["full"], control_rule=CONTROL_RULE["full"],
                                 alias=alias, seed=None, selection=None, members=None,
                                 control_sources=None)
