"""Resolve sealed experiment source identities at their accepted Git revisions."""

from functools import lru_cache
from hashlib import sha256
import subprocess

from cpst_util.release_names import (historical_path, original_alpha_release_sha256,
                                     original_axis_manifest_sha256)


@lru_cache
def git_code_hashes(repo, revision, paths):
    """Hash the exact Git blobs named by a registration."""
    roots = subprocess.check_output(["git", "ls-tree", "--name-only", revision], cwd=repo, text=True).splitlines()
    result = {}
    for path in paths:
        blob_path = historical_path(path) if path.startswith("cpst/") and "cpst" not in roots else path
        contents = subprocess.check_output(["git", "cat-file", "blob", f"{revision}:{blob_path}"], cwd=repo)
        result[path] = sha256(contents).hexdigest()
    return result


def legacy_configuration(config, repo, revision, *, historical_paths=None):
    """Reconstruct the configuration before the shared runtime recorded its identity."""
    original = {name: value for name, value in config.items() if name != "experiment"}
    paths = historical_paths if historical_paths is not None else tuple(
        path for path in config["code_sha256"] if "/util/training/" not in path
    )
    original["code_sha256"] = git_code_hashes(repo, revision, paths)
    return original


def verify_sealed_code(config, repo, revision, *, historical_paths=None):
    if "experiment" in config:
        return
    recorded = config["code_sha256"]
    paths = historical_paths if historical_paths is not None else tuple(sorted(recorded))
    original = git_code_hashes(repo, revision, paths)
    if not recorded or recorded != original:
        raise ValueError(f"sealed registration source differs from Git revision {revision}")


def verify_replay_configuration(recorded, live, *, replay, repo, revision, historical_paths=None):
    if recorded == live:
        return
    if replay and recorded == legacy_configuration(live, repo, revision, historical_paths=historical_paths):
        return
    raise ValueError("registered study, inputs, source or runtime changed")


def historical_gate_configuration(config):
    """Restore the representation of a normalized configuration hashed by a pre-rename gate."""
    prior = {**config, "code_sha256": {historical_path(path): digest
                                      for path, digest in config["code_sha256"].items()}}
    for name in ("schema", "protocol"):
        if name in prior:
            prior[name] = prior[name].replace("cpst", "carps-vae")
    if "axis_family_manifest_sha256" in prior:
        prior["axis_family_manifest_sha256"] = original_axis_manifest_sha256(prior["axis_family_manifest_sha256"])
    if "alpha_release_sha256" in prior:
        prior["alpha_release_sha256"] = original_alpha_release_sha256()
    return prior
