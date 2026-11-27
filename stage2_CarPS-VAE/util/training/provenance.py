"""Resolve sealed experiment source identities at their accepted Git revisions."""

from functools import lru_cache
from hashlib import sha256
import subprocess


@lru_cache
def git_code_hashes(repo, revision, paths):
    """Hash the exact Git blobs named by a registration."""
    return {
        path: sha256(subprocess.check_output(["git", "cat-file", "blob", f"{revision}:{path}"], cwd=repo)).hexdigest()
        for path in paths
    }


def legacy_configuration(config, repo, revision):
    """Reconstruct the configuration before the shared runtime recorded its identity."""
    original = {name: value for name, value in config.items() if name != "experiment"}
    paths = tuple(path for path in config["code_sha256"] if "/util/training/" not in path)
    original["code_sha256"] = git_code_hashes(repo, revision, paths)
    return original


def verify_sealed_code(config, repo, revision):
    if "experiment" in config:
        return
    recorded = config["code_sha256"]
    original = git_code_hashes(repo, revision, tuple(sorted(recorded)))
    if not recorded or recorded != original:
        raise ValueError(f"sealed registration source differs from Git revision {revision}")


def verify_replay_configuration(recorded, live, *, replay, repo, revision):
    if recorded == live:
        return
    if replay and recorded == legacy_configuration(live, repo, revision):
        return
    raise ValueError("registered study, inputs, source or runtime changed")
