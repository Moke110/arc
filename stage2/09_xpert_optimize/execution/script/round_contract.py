"""Round 1 artifact identity and committed validation-promotion gate."""

import hashlib
import json
from pathlib import Path
import subprocess


PROMOTION_SCHEMA = "step09-promotion-v1"


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def checked_payload(root, manifest, relative):
    root = Path(root).absolute()
    path = root / relative
    if (not path.is_relative_to(root) or path.resolve() != path
            or sha256(path) != manifest["files"][relative]):
        raise ValueError(f"artifact path or digest differs from its manifest: {path}")
    return path


def verify_committed(path):
    path = Path(path).absolute()
    repo = next((parent for parent in path.parents if (parent / ".git").exists()), None)
    if repo is None or path.resolve() != path:
        raise ValueError("artifact must have committed dev bytes")
    expected = subprocess.run(
        ["git", "-C", str(repo), "show", f"dev:{path.relative_to(repo).as_posix()}"],
        capture_output=True, check=False,
    )
    if expected.returncode or expected.stdout != path.read_bytes():
        raise ValueError(f"artifact must match committed dev bytes: {path}")
    return expected.stdout


def verify_committed_promotion(path):
    promotion = json.loads(verify_committed(path))
    if (promotion.get("schema"), promotion.get("fold"), promotion.get("seed")) != (PROMOTION_SCHEMA, 0, 17):
        raise ValueError("invalid committed Round 1 promotion")
    for field in ("plan_sha256", "validation_manifest_sha256", "validation_labels_sha256",
                  "validation_metrics_sha256"):
        digest = promotion.get(field, "")
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError(f"invalid committed Round 1 promotion digest: {field}")
    if not isinstance(promotion.get("comparisons"), list) or len(promotion["comparisons"]) != 8:
        raise ValueError("committed Round 1 promotion must account for all eight paired comparisons")
    return promotion
