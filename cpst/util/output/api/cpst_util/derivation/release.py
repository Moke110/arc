"""Candidate-then-replace publication harness for accepted Stage 2 releases.

ADR 0062: a build writes a complete candidate **outside** the accepted
``output/`` tree, verifies it, and only then replaces the accepted release, so a
failed build leaves the accepted release usable. Determinism is byte integrity:
replaying a build from the same pinned inputs must reproduce identical output
digests, separately from any tolerance on stochastic results.

The harness owns only publication mechanics (staging, digest checks, atomic
replacement, manifest path hygiene). Product schemas and science belong to the
product modules that call it (``raw_counts.py``, ``corpus.py``).
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

CHUNK = 8 * 1024 * 1024


class ReleaseError(ValueError):
    """A candidate, manifest or replay violates the release contract."""


def sha256(path):
    """Streaming sha256 of one file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return path


def read_json(path):
    return json.loads(Path(path).read_text())


def output_digests(manifest):
    """Map every declared output path to its sha256."""
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ReleaseError("manifest declares no outputs")
    digests = {}
    for name, record in outputs.items():
        digest = record.get("sha256") if isinstance(record, dict) else record
        if not isinstance(digest, str) or len(digest) != 64:
            raise ReleaseError(f"invalid output digest for {name}")
        digests[name] = digest
    return digests


def check_outputs(root, manifest):
    """Re-hash every declared output under ``root``; fail on any mismatch."""
    root = Path(root)
    checked = {}
    for name, digest in sorted(output_digests(manifest).items()):
        path = root / name
        if not path.is_file():
            raise ReleaseError(f"missing release payload: {name}")
        actual = sha256(path)
        if actual != digest:
            raise ReleaseError(f"release payload digest differs: {name}")
        checked[name] = actual
    return checked


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from _strings(key)
            yield from _strings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _strings(item)


def check_relative_paths(manifest):
    """No release manifest may leak an absolute host path."""
    for text in _strings(manifest):
        if text.startswith("/") or (len(text) > 2 and text[1] == ":" and text[2] in "\\/"):
            raise ReleaseError(f"manifest contains an absolute path: {text}")


def check_manifest(root, manifest):
    """Structural check shared by every product: outputs exist and paths are relative."""
    check_relative_paths(manifest)
    return check_outputs(root, manifest)


def stage_root(output, root=None):
    """Directory for candidate payloads; ADR 0062 keeps build state out of ``output/``."""
    if root is not None:
        return Path(root)
    return Path(output).resolve().parent.parent / "execution" / "checkpoints"


def candidate(output, root=None):
    """Create an empty candidate directory outside ``output/``."""
    stage = stage_root(output, root)
    stage.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f".{Path(output).name}-candidate-", dir=stage))


def replace(candidate_dir, output):
    """Atomically swap a verified candidate into the accepted release path."""
    candidate_dir, output = Path(candidate_dir), Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not output.exists():
        os.rename(candidate_dir, output)
        return output
    previous = output.with_name(f".{output.name}.previous")
    if previous.exists():
        shutil.rmtree(previous)
    os.rename(output, previous)
    try:
        os.rename(candidate_dir, output)
    except BaseException:
        os.rename(previous, output)
        raise
    shutil.rmtree(previous)
    return output


def publish(output, build, verify, root=None, keep_failed=False):
    """Build and verify a candidate outside ``output/``, then replace the release."""
    output = Path(output)
    staging = candidate(output, root)
    try:
        build(staging)
        verify(staging)
        replace(staging, output)
    except BaseException:
        if not keep_failed:
            shutil.rmtree(staging, ignore_errors=True)
        raise
    return output


def replay(output, build, verify, root=None):
    """Rebuild from pinned inputs into a scratch candidate and compare every digest.

    This is the determinism gate: the replayed manifest and every declared output
    digest must equal the accepted release. The accepted release is never written.
    """
    output = Path(output)
    accepted = read_json(output / "manifest.json")
    check_manifest(output, accepted)
    staging = candidate(output, root)
    try:
        build(staging)
        verify(staging)
        replayed = read_json(staging / "manifest.json")
        if replayed != accepted:
            raise ReleaseError("replayed manifest differs from the accepted release")
        check_manifest(staging, accepted)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return {"result": "PASS", "manifest_sha256": sha256(output / "manifest.json"),
            "outputs": sorted(output_digests(accepted))}
