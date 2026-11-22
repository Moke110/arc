#!/usr/bin/env python3
"""Sealed acceptance receipt for the `tsg_view` study (schema `tsg-view-acceptance-v1`).

The builder is a read-only auditor: it hashes the accepted products, parses the capped
attempt telemetry, re-runs the figure renderer from ``figure/data/*.csv`` alone in a fresh
process, proves that the tracked files and named sealed artifacts of Steps 01-06 still match
their pinned baseline, and refuses (fail closed) on a missing artifact, a digest mismatch, a
cap violation, a non-identical re-render or a changed preceding step.

Two phases:

* ``--phase partial`` -- the 45-fit matrix is still running.  Missing *pending* products
  (not-yet-run origins, the Tahoe report pass, the sealed matrix verification) are reported
  as completeness, with the missing origins named.  Corrupt or drifted accepted products
  still fail.
* ``--phase final`` -- the sealed receipt.  Every declared product must be present and
  verified; ``--verify-binary`` is mandatory, so multi-gigabyte payloads are re-hashed
  against their manifest digests; the pre-registration block digest and the ``verify_all.sh``
  transcript are required.

Determinism: the receipt bytes are a pure function of the inputs, except ``generated_at``,
which is injectable (``--timestamp``) so tests can prove byte equality.  All recorded paths
are repo-relative; no absolute ``/tmp`` path can leak into the receipt.

CLI::

    acceptance.py health --root <repo> --output <health-dir> [--verify-all]
    acceptance.py receipt --root <repo> --phase partial|final
                         [--runs DIR] [--truth DIR] [--views DIR] [--freezes DIR]
                         [--inputs DIR] [--report DIR] [--figure DIR] [--log DIR]
                         [--health-dir DIR] [--output FILE]
                         [--baseline FILE] [--refs FILE] [--baseline-receipt FILE]
                         [--preregistration FILE] [--verify-binary] [--timestamp ISO]
                         [--no-figure-replay]
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = "tsg-view-acceptance-v1"
HEALTH_SCHEMA = "tsg-view-health-v1"
INVARIANCE_SCHEMA = "tsg-view-steps01-06-invariance-v1"
PRODUCT_SCHEMA = "tsg-view-product-v1"
FREEZE_SCHEMA = "tsg-view-freeze-v1"
TRUTH_SCHEMA = "tsg-view-truth-v1"
TRUTH_FOLD_SCHEMA = "tsg-view-truth-fold-v1"
RUN_SCHEMA = "tsg-view-run-v1"
RECORD_SCHEMA = "tsg-view-run-record-v1"
SELECTION_SCHEMA = "tsg-view-validation-selection-v1"
PLAN_SCHEMA = "tsg-view-formal-plan-v1"
PLAN_MANIFEST_SCHEMA = "tsg-view-formal-plan-manifest-v1"
CAP_SCHEMA = "tsg-view-formal-cap-v1"
VERIFY_SCHEMA = "tsg-view-formal-verify-v1"
REPORT_SCHEMA = "tsg-view-report-v1"
RENDER_SCHEMA = "tsg-view-render-v1"
SMOKE_SCHEMA = "tsg-view-smoke-v1"

STUDY = "stage2/tsg_view"
SCENARIO = "cold-drug-dose"
SOURCES = ("tenv", "tahoe_bulk")
MAIN_VIEWS = ("full", "merge", "solo", "full-pooled-label")
VIEW_PRODUCTS = {
    "tenv": ("full", "merge", "solo", "solo-deepest", "full-pooled-label"),
    "tahoe_bulk": ("full", "merge", "solo", "full-pooled-label"),
}
VIEW_ARMS = {
    "tenv": ("full", "merge", "solo", "full-pooled-label", "solo-deepest"),
    "tahoe_bulk": ("full", "merge", "solo"),
}
FOLDS = (0, 1, 2)
SEEDS = (42, 43)
SOLO_DEEPEST_SEEDS = (42,)
SOLO_SEED = 20261122
EXPECTED_FITS = 45
MEMORY_MAX = 10 * 1024 ** 3
SWAP_MAX = 8 * 1024 ** 3
CAP_TEXT = ("10G", "8G")
TELEMETRY_SECONDS = 10.0
# the sampler sleeps 10 s between samples; measured jitter reaches 14.66 s, so
# one missed sample (20 s) is the stalled-sampler threshold
TELEMETRY_TOLERANCE = TELEMETRY_SECONDS

FREEZE_REQUIRED = ("config.json", "population.csv", "all_obs.csv", "parents.csv",
                   "genes.csv", "splits.csv", "matched_controls.json", "expression.npy",
                   "control.npy", "morgan.npy", "control_sources.csv")
STEPS = ("01_bases", "02_models_tenv-5k", "03_tsg_drug_emb", "04_tsg_dose",
         "05_tsg_width_dose", "06_tsg_transfer")
BASELINE_REL = (".scratch/tsg-view/research/baseline/"
                "steps01-06-per-path.sha256")
REFS_REL = ".scratch/tsg-view/research/baseline/refs.txt"

VERIFY_ALL_NOTE = {
    "status": "documented-not-run",
    "reason": "heavy offline job (52 GB hashed, ~4 min 33 s CPU, 9.7 GiB page cache); "
              "it must run alone, never beside a fit",
    "reference": "docs/agents/execution-environment.md#known-heavy-jobs",
    "command": ["bash", "script/verify_all.sh"],
}

_BASELINE_LINE = re.compile(r"^([0-9a-f]{64})  (\S+)$")
_RSS = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")
_EXIT = re.compile(r"Exit status: (\d+)")
_UNITS = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4}


class AcceptanceError(RuntimeError):
    """A sealed acceptance invariant did not hold."""


class MissingArtifact(AcceptanceError):
    pass


class DigestMismatch(AcceptanceError):
    pass


class CapViolation(AcceptanceError):
    pass


class InvarianceChanged(AcceptanceError):
    pass


# --------------------------------------------------------------------------- #
# primitives
# --------------------------------------------------------------------------- #
def sha256_file(path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return path


def write_receipt(receipt, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(receipt_bytes(receipt))
    return path


def receipt_bytes(receipt) -> bytes:
    return (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode()


def rel(root, path) -> str:
    root = Path(root).resolve()
    path = Path(path).resolve()
    try:
        return path.relative_to(root).as_posix()
    except ValueError as error:
        raise AcceptanceError(f"path escapes the study root: {path.name}") from error


def file_record(root, path) -> dict:
    path = Path(path)
    return {"sha256": sha256_file(path), "bytes": path.stat().st_size}


def _pin(inventory, group, root, path, record=None):
    path = Path(path)
    if not path.is_file():
        raise MissingArtifact(f"missing artifact: {rel(root, path)}")
    inventory.setdefault(group, {})[rel(root, path)] = record or file_record(root, path)
    return path


def _declared(path, declared, *, verify_binary, inventory, group, root, label):
    """Pin a file whose digest the producing manifest already declares."""
    path = Path(path)
    if not path.is_file():
        raise MissingArtifact(f"missing {label}: {rel(root, path)}")
    digest = declared.get("sha256")
    if not re.fullmatch(r"[0-9a-f]{64}", str(digest)):
        raise AcceptanceError(f"{label} declares no sha256: {rel(root, path)}")
    if verify_binary:
        actual = file_record(root, path)
        if actual["sha256"] != digest:
            raise DigestMismatch(f"digest mismatch: {label} {rel(root, path)}")
        declared_bytes = declared.get("bytes")
        if declared_bytes is not None and actual["bytes"] != declared_bytes:
            raise DigestMismatch(f"digest mismatch: {label} size {rel(root, path)}")
        record = actual
    else:
        record = {"sha256": digest, "bytes": path.stat().st_size}
    inventory.setdefault(group, {})[rel(root, path)] = record
    return path


def _label_for(root, path) -> str:
    try:
        return rel(root, path)
    except AcceptanceError:
        return f"{STUDY}/log/acceptance/health/{Path(path).name}"


def _pin_health(inventory, root, path):
    """Pin a health artifact even when the transcript directory lives outside the root."""
    path = Path(path)
    try:
        key = rel(root, path)
    except AcceptanceError:
        key = f"{STUDY}/log/acceptance/health/{path.name}"
    inventory.setdefault("health", {})[key] = file_record(root, path)


def parse_memory(value) -> int:
    if isinstance(value, int):
        return value
    match = re.fullmatch(r"(\d+)([KMGT])?B?", str(value).strip().upper())
    if not match:
        raise AcceptanceError(f"cannot parse a memory cap: {value!r}")
    return int(match.group(1)) * _UNITS[match.group(2) or ""]


def default_paths(root):
    root = Path(root)
    step = root / STUDY
    return {
        "root": root,
        "step": step,
        "views": step / "output/views",
        "freezes": step / "execution/freezes",
        "truth": step / "output/truth",
        "runs": step / "output/runs",
        "report": step / "output/report",
        "figure": step / "figure",
        "log": step / "log",
        "health": step / "log/acceptance/health",
        "output": step / "log/acceptance/receipt.json",
        "inputs": root / "stage2/01_bases/output/view-5k",
        "baseline": root / BASELINE_REL,
        "refs": root / REFS_REL,
    }


# --------------------------------------------------------------------------- #
# plan / identity
# --------------------------------------------------------------------------- #
def plan_origins():
    """The sealed 45-fit schedule: every tenv fit first, then every Tahoe fit."""
    rows = []
    for source in SOURCES:
        for view in VIEW_ARMS[source]:
            seeds = SOLO_DEEPEST_SEEDS if view == "solo-deepest" else SEEDS
            for fold in FOLDS:
                for seed in seeds:
                    rows.append({"source": source, "view": view, "fold": fold,
                                 "seed": seed, "method": "no-dose",
                                 "scenario": SCENARIO})
    return rows


def origin_key(origin) -> str:
    return (f"{origin['source']}/{origin['view']}/fold-{origin['fold']}"
            f"/seed-{origin['seed']}")


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def git_revs(root) -> dict:
    values = {}
    for name, ref in (("dev", "dev"), ("main", "main"), ("head", "HEAD")):
        completed = subprocess.run(["git", "-C", str(root), "rev-parse", "--verify", ref],
                                   capture_output=True, text=True)
        values[name] = completed.stdout.strip() if completed.returncode == 0 else None
    branch = subprocess.run(["git", "-C", str(root), "rev-parse", "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True)
    values["branch"] = branch.stdout.strip() if branch.returncode == 0 else None
    return values


def _dist_version(name):
    try:
        return importlib.metadata.version(name)
    except Exception:  # pragma: no cover - environment dependent
        return None


def environment_probe(python=None) -> dict:
    python = python or sys.executable
    program = (
        "import json, sys\n"
        "data = {'python': sys.version.split()[0]}\n"
        "try:\n"
        "    import torch\n"
        "    data['torch'] = torch.__version__\n"
        "    data['cuda'] = torch.version.cuda\n"
        "    available = bool(torch.cuda.is_available())\n"
        "    data['device'] = 'cuda' if available else 'cpu'\n"
        "    data['gpu'] = torch.cuda.get_device_name(0) if available else None\n"
        "except Exception as error:\n"
        "    data['torch_error'] = type(error).__name__\n"
        "try:\n"
        "    import numpy\n"
        "    data['numpy'] = numpy.__version__\n"
        "except Exception:\n"
        "    pass\n"
        "print(json.dumps(data))\n"
    )
    try:
        completed = subprocess.run([python, "-I", "-c", program], capture_output=True,
                                   text=True, timeout=180)
        values = json.loads(completed.stdout) if completed.returncode == 0 else {}
    except Exception:  # pragma: no cover - environment dependent
        values = {}
    values.setdefault("python", platform.python_version())
    values.setdefault("device", "cpu")
    driver = None
    try:
        probe = subprocess.run(["nvidia-smi", "--query-gpu=driver_version",
                                "--format=csv,noheader"], capture_output=True, text=True,
                               timeout=30)
        if probe.returncode == 0 and probe.stdout.strip():
            driver = probe.stdout.splitlines()[0].strip()
    except Exception:  # pragma: no cover - environment dependent
        driver = None
    values["driver"] = driver
    values["matplotlib"] = _dist_version("matplotlib")
    return {key: values.get(key) for key in
            ("python", "torch", "cuda", "device", "driver", "numpy", "matplotlib", "gpu")}


# --------------------------------------------------------------------------- #
# telemetry and caps
# --------------------------------------------------------------------------- #
def summarise_phase(phase_dir) -> dict:
    """Summarise one capped driver phase: cap, 10 s samples, peaks, exit status."""
    phase_dir = Path(phase_dir)
    phase = load_json(phase_dir / "phase.json")
    result_path = phase_dir / "result.json"
    result = load_json(result_path) if result_path.is_file() else {}
    time_text = (phase_dir / "time.txt").read_text() if (phase_dir / "time.txt").is_file() else ""
    rss_match = _RSS.search(time_text)
    if rss_match is None:
        raise MissingArtifact(f"missing RSS evidence in phase {phase_dir.name}")
    exit_match = _EXIT.search(time_text)
    time_exit = int(exit_match.group(1)) if exit_match else None
    rows = []
    path = phase_dir / "telemetry.jsonl"
    if path.is_file():
        for line in path.read_text().splitlines():
            if line.strip():
                rows.append(json.loads(line))
    if not rows:
        raise CapViolation(f"cap violation: no 10 s telemetry samples in phase {phase_dir.name}")
    times = [float(row["utc_unix_seconds"]) for row in rows]
    gaps = [right - left for left, right in zip(times, times[1:])]
    if any(gap < 0 for gap in gaps):
        raise CapViolation(f"cap violation: non-monotonic telemetry in phase {phase_dir.name}")
    memory_max = phase.get("memory_max")
    swap_max = phase.get("memory_swap_max")
    return {
        "phase": phase.get("phase"),
        "action": phase.get("action"),
        "started_utc_unix_seconds": phase.get("started_utc_unix_seconds"),
        "threads": phase.get("threads"),
        "memory_max": memory_max,
        "memory_swap_max": swap_max,
        "memory_max_bytes": parse_memory(memory_max),
        "swap_max_bytes": parse_memory(swap_max),
        "samples": len(rows),
        "first_utc_unix_seconds": times[0],
        "last_utc_unix_seconds": times[-1],
        "max_gap_seconds": round(max(gaps), 3) if gaps else 0.0,
        "min_available_ram_kib": min(int(row["mem_available_kib"]) for row in rows),
        "min_free_swap_kib": min(int(row["swap_free_kib"]) for row in rows),
        "peak_gpu_used_mib": max(int(row["gpu_used_mib"]) for row in rows),
        "max_gpu_util_pct": max(int(row.get("gpu_util_pct") or 0) for row in rows),
        "peak_rss_kib": int(rss_match.group(1)),
        "time_exit_status": time_exit,
        "returncode": result.get("returncode"),
        "elapsed_seconds": result.get("elapsed_seconds"),
    }


def assert_caps(summary, *, telemetry_seconds=TELEMETRY_SECONDS) -> dict:
    if str(summary.get("memory_max")) != CAP_TEXT[0] or \
            summary.get("memory_max_bytes") != MEMORY_MAX:
        raise CapViolation(
            f"cap violation: MemoryMax {summary.get('memory_max')!r} != {CAP_TEXT[0]}")
    if str(summary.get("memory_swap_max")) != CAP_TEXT[1] or \
            summary.get("swap_max_bytes") != SWAP_MAX:
        raise CapViolation(
            f"cap violation: MemorySwapMax {summary.get('memory_swap_max')!r} != "
            f"{CAP_TEXT[1]}")
    peak = summary.get("peak_rss_kib")
    if peak is None:
        raise CapViolation("cap violation: no peak RSS measurement")
    if int(peak) * 1024 > MEMORY_MAX:
        raise CapViolation(f"cap violation: peak RSS {peak} KiB exceeds MemoryMax")
    if int(summary.get("samples") or 0) < 1:
        raise CapViolation("cap violation: no telemetry samples")
    if float(summary.get("max_gap_seconds") or 0.0) > telemetry_seconds + TELEMETRY_TOLERANCE:
        raise CapViolation(
            f"cap violation: telemetry cadence gap "
            f"{summary.get('max_gap_seconds')} s exceeds {telemetry_seconds:g} s")
    if summary.get("returncode") not in (0, None) or \
            summary.get("time_exit_status") not in (0, None):
        raise AcceptanceError(
            f"phase exit status {summary.get('returncode')}/"
            f"{summary.get('time_exit_status')} on {summary.get('action')!r}")
    return summary


def verify_origin_resources(root, resources_root, origin, *, phase, inventory):
    folder = (Path(resources_root) / origin["source"] / origin["view"]
              / f"fold-{origin['fold']}" / f"seed-{origin['seed']}")
    if not folder.is_dir():
        if phase == "final":
            raise MissingArtifact(f"missing telemetry for {origin_key(origin)}")
        return None
    attempts = sorted(path for path in folder.iterdir()
                      if path.is_dir() and (path / "phase.json").is_file())
    if not attempts:
        raise MissingArtifact(f"missing telemetry phases for {origin_key(origin)}")
    summaries = {}
    for attempt in attempts:
        summary = assert_caps(summarise_phase(attempt))
        for name in ("phase.json", "telemetry.jsonl", "time.txt", "result.json"):
            _pin(inventory, "resources", root, attempt / name)
        summaries[attempt.name] = {
            "phase": summary["phase"], "action": summary["action"],
            "samples": summary["samples"], "max_gap_seconds": summary["max_gap_seconds"],
            "peak_rss_kib": summary["peak_rss_kib"],
            "peak_gpu_used_mib": summary["peak_gpu_used_mib"],
            "min_available_ram_kib": summary["min_available_ram_kib"],
            "returncode": summary["returncode"],
            "elapsed_seconds": summary["elapsed_seconds"],
            "memory_max_bytes": summary["memory_max_bytes"],
            "swap_max_bytes": summary["swap_max_bytes"],
        }
    actions = [value["action"] for value in summaries.values()]
    if "train" not in actions or "evaluate" not in actions or actions.count("replay") < 2:
        raise AcceptanceError(
            f"replay evidence missing for {origin_key(origin)}: "
            f"expected train, evaluate and two fresh-process replays, got {actions}")
    return {
        "origin": origin_key(origin),
        "phases": summaries,
        "telemetry_samples": sum(value["samples"] for value in summaries.values()),
        "max_sample_gap_seconds": max(value["max_gap_seconds"]
                                      for value in summaries.values()),
        "peak_rss_kib": max(value["peak_rss_kib"] for value in summaries.values()),
        "peak_gpu_used_mib": max(value["peak_gpu_used_mib"] for value in summaries.values()),
        "min_available_ram_kib": min(value["min_available_ram_kib"]
                                     for value in summaries.values()),
        "memory_max_bytes": MEMORY_MAX,
        "swap_max_bytes": SWAP_MAX,
    }


# --------------------------------------------------------------------------- #
# view products / freezes / truth
# --------------------------------------------------------------------------- #
def verify_views(root, views_root, *, phase, verify_binary, inventory):
    products = {}
    for source in SOURCES:
        for view in VIEW_PRODUCTS[source]:
            path = Path(views_root) / view / f"{source}.json"
            if not path.is_file():
                raise MissingArtifact(f"missing view manifest: {rel(root, path)}")
            manifest = load_json(path)
            if manifest.get("schema") != PRODUCT_SCHEMA:
                raise AcceptanceError(
                    f"view manifest schema {manifest.get('schema')!r} != {PRODUCT_SCHEMA}")
            if manifest.get("source") != source or manifest.get("view") != view:
                raise AcceptanceError(
                    f"view manifest identity mismatch: {rel(root, path)}")
            _pin(inventory, "views", root, path)
            axis = manifest["axis"]
            axis_path = Path(root) / axis["file"]
            if sha256_file(axis_path) != axis["sha256"]:
                raise DigestMismatch(f"digest mismatch: gene axis {rel(root, axis_path)}")
            _pin(inventory, "views", root, axis_path)
            sources = manifest["source_manifests"]
            release = Path(root) / sources["accepted_view"]
            release_digest = sha256_file(release)
            if release_digest != sources["accepted_view_sha256"]:
                raise DigestMismatch(f"digest mismatch: accepted view release "
                                     f"{rel(root, release)}")
            _pin(inventory, "views", root, release)
            payload = manifest.get("payload")
            alias = manifest.get("alias")
            payload_record = None
            if payload:
                payload_path = path.parent / payload["file"]
                _declared(payload_path, payload, verify_binary=verify_binary,
                          inventory=inventory, group="view_payloads", root=root,
                          label="view payload")
                payload_record = {"file": rel(root, payload_path),
                                  "sha256": payload["sha256"],
                                  "bytes": payload.get("bytes"),
                                  "kind": payload.get("kind", "payload"),
                                  "rows": payload.get("rows"), "cols": payload.get("cols")}
            elif not alias or not alias.get("manifest_sha256"):
                raise AcceptanceError(
                    f"view product has neither payload nor alias: {rel(root, path)}")
            if alias and alias.get("manifest_sha256"):
                if alias["manifest_sha256"] != release_digest or \
                        sources["accepted_view_sha256"] != release_digest:
                    raise DigestMismatch(f"digest mismatch: view alias {rel(root, path)}")
            companions = {}
            for field in ("profile_ids", "members", "control_sources"):
                declared = manifest.get(field)
                if not declared:
                    continue
                companion = path.parent / declared["file"]
                if not companion.is_file():
                    raise MissingArtifact(f"missing view companion: {rel(root, companion)}")
                if sha256_file(companion) != declared["sha256"]:
                    raise DigestMismatch(
                        f"digest mismatch: view companion {rel(root, companion)}")
                _pin(inventory, "view_companions", root, companion)
                companions[field] = declared["sha256"]
            key = f"{source}-{view}"
            products[key] = {
                "product": manifest.get("product"), "source": source, "view": view,
                "manifest": rel(root, path), "manifest_sha256": sha256_file(path),
                "rows": manifest.get("rows"), "conditions": manifest.get("conditions"),
                "seed": manifest.get("seed"), "selection": manifest.get("selection"),
                "payload": payload_record, "mode": "alias" if payload is None else "payload",
                "companions": companions,
            }
    expected = sum(len(views) for views in VIEW_PRODUCTS.values())
    if len(products) != expected:
        raise MissingArtifact(f"missing view products: {expected - len(products)}")
    return products


def verify_freezes(root, freezes_root, view_products, *, phase, verify_binary, inventory):
    products = {}
    for source in SOURCES:
        for view in VIEW_PRODUCTS[source]:
            folder = Path(freezes_root) / view / source
            path = folder / "manifest.json"
            if not path.is_file():
                raise MissingArtifact(f"missing freeze manifest: {rel(root, path)}")
            manifest = load_json(path)
            if manifest.get("schema") != FREEZE_SCHEMA:
                raise AcceptanceError(
                    f"freeze manifest schema {manifest.get('schema')!r} != {FREEZE_SCHEMA}")
            _pin(inventory, "freezes", root, path)
            outputs = manifest.get("outputs") or {}
            missing = [name for name in FREEZE_REQUIRED if name not in outputs]
            if missing:
                raise MissingArtifact(
                    f"missing freeze outputs in {rel(root, path)}: {', '.join(missing)}")
            for name, digest in sorted(outputs.items()):
                _declared(folder / name, {"sha256": digest}, verify_binary=verify_binary,
                          inventory=inventory, group="freezes", root=root,
                          label="freeze output")
            product = manifest.get("view_product") or {}
            view_sha = view_products[f"{source}-{view}"]["manifest_sha256"]
            if product.get("sha256") != view_sha:
                raise DigestMismatch(
                    f"digest mismatch: freeze view product {rel(root, path)}")
            sealed = manifest.get("sealed") or {}
            sealed_files = {}
            for name, key in (("population.csv", "population_sha256"),
                              ("splits.csv", "splits_sha256"),
                              ("parents.csv", "parents_sha256"),
                              ("morgan.npy", "morgan_sha256")):
                sealed_path = Path(root) / sealed["file"] / name
                if not sealed_path.is_file():
                    raise MissingArtifact(f"missing sealed input: {rel(root, sealed_path)}")
                if sha256_file(sealed_path) != sealed[key]:
                    raise DigestMismatch(
                        f"digest mismatch: sealed input {rel(root, sealed_path)}")
                sealed_files[name] = sealed[key]
            products[f"{source}-{view}"] = {
                "manifest": rel(root, path), "manifest_sha256": sha256_file(path),
                "delta_kind": manifest.get("delta_kind"), "outputs": len(outputs),
                "view_product_mode": product.get("mode"), "sealed": sealed_files,
            }
    expected = sum(len(views) for views in VIEW_PRODUCTS.values())
    if len(products) != expected:
        raise MissingArtifact(f"missing freezes: {expected - len(products)}")
    return products


def verify_truth(root, truth_root, *, phase, verify_binary, inventory):
    products = {}
    for source in SOURCES:
        folder = Path(truth_root) / source
        path = folder / "manifest.json"
        if not path.is_file():
            raise MissingArtifact(f"missing truth manifest: {rel(root, path)}")
        manifest = load_json(path)
        if manifest.get("schema") != TRUTH_SCHEMA:
            raise AcceptanceError(
                f"truth manifest schema {manifest.get('schema')!r} != {TRUTH_SCHEMA}")
        if manifest.get("source") != source:
            raise AcceptanceError(f"truth manifest identity mismatch: {rel(root, path)}")
        _pin(inventory, "truth", root, path)
        for key in ("diagnostics", "noise_floor"):
            declared = manifest[key]
            target = folder / declared["file"]
            matched = _declared(target, declared, verify_binary=verify_binary,
                                inventory=inventory, group="truth", root=root,
                                label=f"truth {key}")
            if sha256_file(target) != declared["sha256"]:
                raise DigestMismatch(f"digest mismatch: truth {key} "
                                     f"{rel(root, target)}")
            del matched
        folds = {}
        for partition in ("primary", "condition_mean"):
            for fold, entry in sorted(manifest[partition].items()):
                fold_dir = folder / partition / fold
                fold_manifest = fold_dir / "manifest.json"
                if not fold_manifest.is_file():
                    raise MissingArtifact(f"missing truth {partition} manifest: "
                                          f"{rel(root, fold_manifest)}")
                if sha256_file(fold_manifest) != entry["manifest"]["sha256"]:
                    raise DigestMismatch(f"digest mismatch: truth {partition} manifest "
                                         f"{rel(root, fold_manifest)}")
                _pin(inventory, "truth", root, fold_manifest)
                for name, declared in sorted(entry["outputs"].items()):
                    _declared(fold_dir / name, declared, verify_binary=verify_binary,
                              inventory=inventory, group="truth", root=root,
                              label=f"truth {partition} output")
                folds[f"{partition}/{fold}"] = {
                    "manifest_sha256": entry["manifest"]["sha256"],
                    "rows": entry.get("rows"), "conditions": entry.get("conditions")}
        products[source] = {
            "manifest": rel(root, path), "manifest_sha256": sha256_file(path),
            "folds": folds,
            "diagnostics_sha256": manifest["diagnostics"]["sha256"],
            "noise_floor_sha256": manifest["noise_floor"]["sha256"],
            "label": manifest["diagnostics"].get("label"),
        }
    return products


# --------------------------------------------------------------------------- #
# runs, replay, capped driver resources
# --------------------------------------------------------------------------- #
def _walk_smoke_claims(value, trail=(), context=""):
    if isinstance(value, dict):
        scalars = {key: child for key, child in value.items()
                   if not isinstance(child, (dict, list))}
        context = (context + " " + json.dumps(scalars, sort_keys=True)).strip()
        if "phases" in value and ("replay" in value or "exact_replay" in value):
            yield trail, value, context
        for key, child in value.items():
            yield from _walk_smoke_claims(child, trail + (str(key),), context)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _walk_smoke_claims(child, trail + (str(index),), context)


def _smoke_claim(record, label, source=None):
    phases = record.get("phases") or []
    replays = [item for item in phases if item.get("phase") == "replay"]
    if len(replays) < 2 or any(item.get("returncode") != 0 for item in phases):
        raise AcceptanceError(f"replay evidence missing: {label}")
    claim = str(record.get("exact_replay") or record.get("replay") or "")
    if "process" not in claim.lower():
        raise AcceptanceError(f"replay evidence missing from {label}: not a fresh-process "
                              f"claim")
    return {
        "claim": label,
        "evidence": rel(Path(record["_root"]), Path(record["_path"])),
        "sha256": record.get("_sha256"), "source": source or record.get("source"),
        "view": record.get("view"), "fold": record.get("fold"), "seed": record.get("seed"),
        "device": record.get("device"), "exact_replay": claim,
        "replay_phases": len(replays), "rows": record.get("rows"),
        "epochs": record.get("epochs") or record.get("epochs_run"),
    }


def verify_smoke(root, step_dir, *, phase, inventory):
    step_dir = Path(step_dir)
    evidence = []
    receipts = sorted((step_dir / "execution/checkpoints").glob("merger-*/smoke-*.json"))
    for path in receipts:
        value = load_json(path)
        if value.get("schema") != SMOKE_SCHEMA:
            raise AcceptanceError(f"smoke receipt schema mismatch: {rel(root, path)}")
        if value.get("status") != "complete":
            raise AcceptanceError(f"smoke run not complete: {rel(root, path)}")
        value.update({"_root": str(root), "_path": str(path),
                      "_sha256": sha256_file(path)})
        _pin(inventory, "smoke", root, path)
        evidence.append(_smoke_claim(value, rel(root, path)))
    for path in sorted((step_dir / "log").glob("build-0*.json")):
        value = load_json(path)
        claims = list(_walk_smoke_claims(value.get("gates") or {}))
        if not claims:
            continue
        _pin(inventory, "smoke", root, path)
        for trail, record, context in claims:
            label = f"{rel(root, path)}:{'.'.join(trail)}"
            record = dict(record)
            record.update({"_root": str(root), "_path": str(path),
                           "_sha256": sha256_file(path)})
            entry = _smoke_claim(record, label)
            if entry["source"] is None:
                haystack = (" ".join(trail) + " " + context + " "
                            + json.dumps(record, sort_keys=True))
                entry["source"] = next((name for name in SOURCES
                                        if name in haystack), None)
            evidence.append(entry)
    if not evidence:
        raise MissingArtifact("missing smoke replay evidence")
    if phase == "final" and not set(SOURCES) <= {item.get("source") for item in evidence}:
        raise MissingArtifact("smoke replay evidence missing for a source")
    return evidence


def verify_runs(root, runs_root, origin_products, freeze_products, truth_products, *,
                phase, verify_binary, inventory):
    runs_root = Path(runs_root)
    plan_path = runs_root / "formal_plan.json"
    plan_manifest_path = runs_root / "formal_plan_manifest.json"
    cap_path = runs_root / "formal_cap.json"
    records_path = runs_root / "run_records.json"
    for path, label in ((plan_path, "formal plan"), (plan_manifest_path, "plan manifest"),
                        (cap_path, "formal cap"), (records_path, "run records")):
        if not path.is_file():
            raise MissingArtifact(f"missing {label}: {rel(root, path)}")
        _pin(inventory, "runs", root, path)
    plan = load_json(plan_path)
    plan_sha = sha256_file(plan_path)
    if plan.get("schema") != PLAN_SCHEMA:
        raise AcceptanceError(f"formal plan schema mismatch: {rel(root, plan_path)}")
    def signature(origin):
        return (origin.get("source"), origin.get("view"), origin.get("fold"),
                origin.get("seed"), origin.get("method"), origin.get("scenario"))
    declared = [signature(origin) for origin in plan.get("origins") or []]
    expected = [signature(origin) for origin in plan_origins()]
    if len(declared) != EXPECTED_FITS or set(declared) != set(expected):
        raise AcceptanceError("plan digest mismatch: declared origins differ from the "
                              "45-fit schedule")
    sources_in_order = [origin.get("source") for origin in plan.get("origins") or []]
    if sources_in_order != sorted(sources_in_order, key=lambda name: SOURCES.index(name)):
        raise AcceptanceError("plan digest mismatch: the schedule does not run every tenv "
                              "fit before the Tahoe fits")
    plan_manifest = load_json(plan_manifest_path)
    if plan_manifest.get("schema") != PLAN_MANIFEST_SCHEMA or \
            plan_manifest.get("plan_sha256") != plan_sha or \
            plan_manifest.get("origins") != EXPECTED_FITS:
        raise AcceptanceError(f"plan digest mismatch: "
                              f"{rel(root, plan_manifest_path)}")
    cap = load_json(cap_path)
    if cap.get("schema") != CAP_SCHEMA or cap.get("memory_max") != CAP_TEXT[0] or \
            cap.get("memory_swap_max") != CAP_TEXT[1]:
        raise CapViolation(f"cap violation: formal cap {rel(root, cap_path)}")
    scope = cap.get("scope") or {}
    if parse_memory(scope.get("MemoryMax")) != MEMORY_MAX or \
            parse_memory(scope.get("MemorySwapMax")) != SWAP_MAX:
        raise CapViolation(f"cap violation: scope caps in {rel(root, cap_path)}")
    if not re.fullmatch(r"[0-9a-f]{32}", str(scope.get("InvocationID"))):
        raise CapViolation(f"cap violation: scope InvocationID in {rel(root, cap_path)}")
    present_plan, pending_origins = [], []
    for origin in plan_origins():
        folder = (runs_root / "runs" / origin["source"] / origin["view"]
                  / f"fold-{origin['fold']}" / f"seed-{origin['seed']}")
        manifest_path = folder / "manifest.json"
        terminal_path = folder / "terminal.json"
        record_path = (runs_root / "records" / origin["source"] / origin["view"]
                       / f"fold-{origin['fold']}" / f"seed-{origin['seed']}.json")
        if not manifest_path.is_file():
            pending_origins.append(origin)
            continue
        if not terminal_path.is_file() or not record_path.is_file():
            if phase == "final":
                raise MissingArtifact(f"missing run artifact: {origin_key(origin)}")
            pending_origins.append(origin)
            continue
        if (load_json(terminal_path).get("status") != "complete"):
            if phase == "final":
                raise AcceptanceError(f"run not complete: {origin_key(origin)}")
            pending_origins.append(origin)
            continue
        present_plan.append(origin)
    needed = {(origin["source"], origin["fold"]) for origin in present_plan}
    if phase == "final":
        needed = {(source, fold) for source in SOURCES for fold in FOLDS}
    selections = {}
    for source, fold in sorted(needed):
        path = runs_root / "selection" / source / f"fold-{fold}.json"
        if not path.is_file():
            raise MissingArtifact(f"missing fold choice receipt: {rel(root, path)}")
        value = load_json(path)
        if value.get("schema") != SELECTION_SCHEMA or \
                value.get("plan_sha256") != plan_sha:
            raise AcceptanceError(f"selection receipt mismatch: {rel(root, path)}")
        _pin(inventory, "runs", root, path)
        selections[f"{source}/fold-{fold}"] = sha256_file(path)
    present, resources, replay = [], {}, {}
    for origin in present_plan:
        key = origin_key(origin)
        folder = (runs_root / "runs" / origin["source"] / origin["view"]
                  / f"fold-{origin['fold']}" / f"seed-{origin['seed']}")
        manifest_path = folder / "manifest.json"
        manifest = load_json(manifest_path)
        if manifest.get("schema") != RUN_SCHEMA:
            raise AcceptanceError(f"run manifest schema mismatch: "
                                  f"{rel(root, manifest_path)}")
        if manifest.get("origin") != origin:
            raise AcceptanceError(f"run identity mismatch: {rel(root, manifest_path)}")
        if manifest.get("plan") != plan_sha:
            raise AcceptanceError(f"plan digest mismatch: {rel(root, manifest_path)}")
        freeze = manifest.get("freeze") or {}
        product = origin_products[f"{origin['source']}-{origin['view']}"]
        if freeze.get("view_product_sha256") != product["manifest_sha256"]:
            raise DigestMismatch(f"digest mismatch: run view product "
                                 f"{rel(root, manifest_path)}")
        freeze = manifest.get("freeze") or {}
        if freeze.get("manifest_sha256") != \
                freeze_products[f"{origin['source']}-{origin['view']}"]["manifest_sha256"]:
            raise DigestMismatch(f"digest mismatch: run freeze "
                                 f"{rel(root, manifest_path)}")
        truth = manifest.get("truth") or {}
        primary = truth.get("primary") or {}
        if primary.get("manifest_sha256") != \
                truth_products[origin["source"]]["folds"][
                    f"primary/fold-{origin['fold']}"]["manifest_sha256"]:
            raise DigestMismatch(f"digest mismatch: run truth "
                                 f"{rel(root, manifest_path)}")
        _pin(inventory, "runs", root, manifest_path)
        for name, digest in sorted((manifest.get("outputs") or {}).items()):
            _declared(folder / name, {"sha256": digest}, verify_binary=verify_binary,
                      inventory=inventory, group="runs", root=root, label="run output")
        terminal_path = folder / "terminal.json"
        if not terminal_path.is_file():
            raise MissingArtifact(f"missing run terminal: {rel(root, terminal_path)}")
        terminal = load_json(terminal_path)
        if terminal.get("status") != "complete":
            raise AcceptanceError(f"run not complete: {rel(root, terminal_path)}")
        claim = str(terminal.get("replay") or "")
        if "fresh-process" not in claim or "exact" not in claim.lower():
            raise AcceptanceError(f"replay evidence missing from "
                                  f"{rel(root, terminal_path)}")
        tokens = {}
        for token in ("validation_replay_token", "test_replay_token",
                      "replay_process_token", "training_process_token",
                      "test_evaluation_process_token"):
            value = terminal.get(token)
            if not re.fullmatch(r"[0-9a-f]{16,64}", str(value)):
                raise AcceptanceError(f"replay evidence missing: {token} in "
                                      f"{rel(root, terminal_path)}")
            tokens[token] = value
        record_path = (runs_root / "records" / origin["source"] / origin["view"]
                       / f"fold-{origin['fold']}" / f"seed-{origin['seed']}.json")
        if not record_path.is_file():
            raise MissingArtifact(f"missing run record: {rel(root, record_path)}")
        recorded = load_json(record_path)
        if recorded.get("schema") != RECORD_SCHEMA or \
                recorded.get("plan_sha256") != plan_sha:
            raise AcceptanceError(f"run record mismatch: {rel(root, record_path)}")
        expected_selection = selections[f"{origin['source']}/fold-{origin['fold']}"]
        if recorded.get("selection_sha256") != expected_selection:
            raise AcceptanceError(f"selection digest mismatch: {rel(root, record_path)}")
        if (recorded.get("artifacts") or {}) != (manifest.get("outputs") or {}):
            raise DigestMismatch(f"digest mismatch: run record artifacts "
                                 f"{rel(root, record_path)}")
        _pin(inventory, "runs", root, record_path)
        resource = verify_origin_resources(root, runs_root / "resources", origin,
                                           phase=phase, inventory=inventory)
        if resource is None:
            resource = {"origin": key, "status": "pending"}
        resources[key] = resource
        replay[key] = {
            "origin": key, "terminal_sha256": sha256_file(terminal_path),
            "validation_replay": "exact", "test_replay": "exact",
            **tokens,
        }
        present.append({
            "origin": origin, "manifest": rel(root, manifest_path),
            "manifest_sha256": sha256_file(manifest_path),
            "record_sha256": sha256_file(record_path),
            "selection_sha256": expected_selection,
            "device": terminal.get("device"),
            "epochs": terminal.get("epochs"), "best_epoch": terminal.get("best_epoch"),
            "validation_mse": terminal.get("validation_mse"),
            "test_samples": terminal.get("test_samples"),
            "peak_rss_kib": recorded.get("peak_rss_kib"),
            "peak_gpu_used_mib": recorded.get("peak_gpu_used_mib"),
            "runtime_seconds": recorded.get("runtime_seconds"),
            "metrics": recorded.get("metrics"),
        })
    records = load_json(records_path)
    if records.get("schema") != "tsg-view-formal-run-v1" or \
            records.get("plan_sha256") != plan_sha or \
            records.get("selected") != EXPECTED_FITS:
        raise AcceptanceError(f"run records mismatch: {rel(root, records_path)}")
    if phase == "final" and (len(present) != EXPECTED_FITS or pending_origins):
        raise AcceptanceError(
            f"incomplete matrix: {len(present)}/{EXPECTED_FITS} origins present, "
            f"{len(pending_origins)} pending")
    verify_path = runs_root / "formal_verify.json"
    verify = None
    if verify_path.is_file():
        value = load_json(verify_path)
        if value.get("schema") != VERIFY_SCHEMA:
            raise AcceptanceError(f"formal verify schema mismatch: "
                                  f"{rel(root, verify_path)}")
        if value.get("status") != "PASS":
            raise AcceptanceError(f"formal verify status {value.get('status')!r} "
                                  f"in {rel(root, verify_path)}")
        if value.get("origins") != EXPECTED_FITS or \
                value.get("plan_sha256") != plan_sha:
            raise AcceptanceError(f"formal verify mismatch: {rel(root, verify_path)}")
        by_origin = {}
        for item in value.get("runs") or []:
            by_origin[json.dumps(item["origin"], sort_keys=True)] = \
                item.get("run_manifest_sha256")
        for entry in present:
            expected = by_origin.get(json.dumps(entry["origin"], sort_keys=True))
            if expected != entry["manifest_sha256"]:
                raise AcceptanceError(f"formal verify run digest mismatch for "
                                      f"{origin_key(entry['origin'])}")
        _pin(inventory, "runs", root, verify_path)
        verify = {"path": rel(root, verify_path), "sha256": sha256_file(verify_path),
                  "status": value.get("status"), "origins": value.get("origins")}
    elif phase == "final":
        raise MissingArtifact(f"missing formal verification report: "
                              f"{rel(root, verify_path)}")
    return {
        "present": present, "pending_origins": pending_origins,
        "resources": resources, "replay": replay,
        "plan_sha256": plan_sha,
        "plan_manifest_sha256": sha256_file(plan_manifest_path),
        "formal_cap": {"sha256": sha256_file(cap_path), "scope": cap.get("scope"),
                       "memory_max_bytes": MEMORY_MAX, "swap_max_bytes": SWAP_MAX},
        "run_records": {"sha256": sha256_file(records_path),
                        "complete": records.get("complete"),
                        "runs": records.get("runs"),
                        "selected": records.get("selected")},
        "formal_verify": verify,
    }


# --------------------------------------------------------------------------- #
# report tables and CSV-only figure reproduction
# --------------------------------------------------------------------------- #
def figure_replay(root, figure_dir, *, python=None, timeout=300) -> dict:
    """Re-render every PNG from ``figure/data/*.csv`` alone and require identical bytes."""
    figure_dir = Path(figure_dir)
    status_path = figure_dir / "render_status.json"
    if not status_path.is_file():
        raise MissingArtifact(f"missing render status: {rel(root, status_path)}")
    status = load_json(status_path)
    if status.get("schema") != RENDER_SCHEMA:
        raise AcceptanceError(f"render status schema mismatch: {rel(root, status_path)}")
    renderer = Path(root) / STUDY / "execution/render.py"
    if not renderer.is_file():
        raise MissingArtifact(f"missing figure renderer: {rel(root, renderer)}")
    data_files = sorted((figure_dir / "data").glob("*.csv"))
    if not data_files:
        raise MissingArtifact("missing figure/data CSV tables")
    python = python or sys.executable
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        (work / "data").mkdir()
        (work / "script").mkdir()
        for path in data_files:
            shutil.copyfile(path, work / "data" / path.name)
        shutil.copyfile(renderer, work / "script" / "render.py")
        environment = dict(os.environ)
        environment.update({"MPLCONFIGDIR": str(work / "mpl"),
                            "PYTHONHASHSEED": "0",
                            "PYTHONDONTWRITEBYTECODE": "1",
                            "HOME": str(work / "home")})
        environment.pop("PYTHONPATH", None)
        (work / "home").mkdir(exist_ok=True)
        completed = subprocess.run(
            [python, "-I", "script/render.py", "--data", "data", "--output", "out",
             "--suite", str(status.get("suite") or "all")],
            cwd=work, env=environment, capture_output=True, text=True, timeout=timeout)
        if completed.returncode != 0:
            raise AcceptanceError(f"figure re-render failed: "
                                  f"{completed.stderr.strip().splitlines()[-1] if completed.stderr.strip() else completed.returncode}")
        produced = {path.name: sha256_file(path)
                    for path in sorted((work / "out").glob("*.png"))}
        if not produced:
            raise AcceptanceError("figure re-render differs: no PNG was produced")
        for name, digest in status.get("figures", {}).items():
            if produced.get(name) != digest:
                raise AcceptanceError(f"figure re-render differs for {name}")
        if set(produced) - set(status.get("figures", {})):
            raise AcceptanceError("figure re-render differs: unexpected PNG")
        replayed = work / "out" / "render_status.json"
        if not replayed.is_file() or replayed.read_bytes() != status_path.read_bytes():
            raise AcceptanceError("figure re-render differs: render_status.json bytes")
    return {
        "suite": status.get("suite"),
        "figures": produced,
        "byte_identical": True,
        "render_status_sha256": sha256_file(status_path),
        "renderer_sha256": sha256_file(renderer),
        "matplotlib": status.get("matplotlib"),
        "dpi": status.get("dpi"),
    }


def verify_report(root, report_root, figure_root, *, phase, inventory, python=None,
                  replay_figures=True):
    report_root = Path(report_root)
    figure_root = Path(figure_root)
    manifest_path = report_root / "report_manifest.json"
    if not manifest_path.is_file():
        if phase == "final":
            raise MissingArtifact(f"missing report manifest: {rel(root, manifest_path)}")
        return {"status": "pending", "manifest": rel(root, manifest_path),
                "tables": {}, "figure_data": {}, "csv_only_replay": {
                    "status": "pending"}, "figures": {}}
    manifest = load_json(manifest_path)
    if manifest.get("schema") != REPORT_SCHEMA:
        raise AcceptanceError(f"report manifest schema mismatch: "
                              f"{rel(root, manifest_path)}")
    _pin(inventory, "report", root, manifest_path)
    tables = {}
    for name, digest in sorted(manifest["tables"].items()):
        path = report_root / name
        if not path.is_file():
            raise MissingArtifact(f"missing report table: {rel(root, path)}")
        if sha256_file(path) != digest:
            raise DigestMismatch(f"digest mismatch: report table {rel(root, path)}")
        _pin(inventory, "report_tables", root, path)
        tables[name] = digest
    on_disk = {path.name for path in report_root.glob("*.csv")}
    if on_disk - set(tables):
        raise AcceptanceError(f"unlisted report table: "
                              f"{', '.join(sorted(on_disk - set(tables)))}")
    figure_data = {}
    for key, digest in sorted(manifest["figure_data"].items()):
        path = figure_root / "data" / Path(key).name
        if not path.is_file():
            raise MissingArtifact(f"missing figure data table: {rel(root, path)}")
        if sha256_file(path) != digest:
            raise DigestMismatch(f"figure data digest mismatch: {rel(root, path)}")
        _pin(inventory, "figure_data", root, path)
        figure_data[path.name] = digest
    if phase == "final":
        if manifest.get("complete_origins") != EXPECTED_FITS or \
                manifest.get("pending_origins") != 0:
            raise AcceptanceError(
                f"incomplete report: {manifest.get('complete_origins')}/"
                f"{EXPECTED_FITS} origins, {manifest.get('pending_origins')} pending")
        if sorted(manifest.get("sources") or []) != sorted(SOURCES):
            raise MissingArtifact("report is missing a source suite")
    status_path = figure_root / "render_status.json"
    figures = {}
    csv_only = {"status": "pending"}
    if status_path.is_file():
        status = load_json(status_path)
        if status.get("schema") != RENDER_SCHEMA:
            raise AcceptanceError(f"render status schema mismatch: "
                                  f"{rel(root, status_path)}")
        _pin(inventory, "figures", root, status_path)
        if sorted(status.get("csv_inputs")) != sorted(figure_data):
            raise AcceptanceError(f"figure inputs differ from the report figure data: "
                                  f"{rel(root, status_path)}")
        for name, digest in sorted(status["csv_inputs"].items()):
            if figure_data.get(name) != digest:
                raise DigestMismatch(f"figure data digest mismatch: {name}")
        for name, digest in sorted(status.get("figures", {}).items()):
            path = figure_root / name
            if not path.is_file():
                raise MissingArtifact(f"missing figure: {rel(root, path)}")
            if sha256_file(path) != digest:
                raise DigestMismatch(f"figure digest mismatch: {rel(root, path)}")
            _pin(inventory, "figures", root, path)
            figures[name] = digest
        on_disk = {path.name for path in figure_root.glob("*.png")}
        if on_disk - set(status.get("figures", {})):
            raise AcceptanceError(f"unlisted figure: "
                                  f"{', '.join(sorted(on_disk - set(figures)))}")
        if phase == "final":
            required = {f"learning_curves_{source}.png" for source in SOURCES}
            required |= {f"metrics_{source}.png" for source in SOURCES}
            required |= {f"metrics_per_dose_{source}.png" for source in SOURCES}
            required |= {"views_effect.png"}
            missing = required - set(figures)
            if missing:
                raise MissingArtifact(f"missing figure: {', '.join(sorted(missing))}")
        if replay_figures:
            csv_only = figure_replay(root, figure_root, python=python)
    elif phase == "final":
        raise MissingArtifact(f"missing render status: {rel(root, status_path)}")
    return {"status": "verified", "manifest": rel(root, manifest_path),
            "manifest_sha256": sha256_file(manifest_path), "tables": tables,
            "figure_data": figure_data, "figures": figures,
            "sources": manifest.get("sources"),
            "complete_origins": manifest.get("complete_origins"),
            "pending_origins": manifest.get("pending_origins"),
            "csv_only_replay": csv_only}


# --------------------------------------------------------------------------- #
# health checks
# --------------------------------------------------------------------------- #
def run_health(root, output_dir, *, scripts=("check_paths", "check_tracked"),
               verify_all=False, timeout=900, commit=None) -> dict:
    root = Path(root)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    names = tuple(scripts)
    if verify_all:
        names = names + ("verify_all",)
    checks = {}
    for name in names:
        script = root / "script" / f"{name}.sh"
        if not script.is_file():
            raise MissingArtifact(f"health script missing: script/{name}.sh")
        started = time.monotonic()
        completed = subprocess.run(["bash", str(script)], cwd=root, capture_output=True,
                                   text=True, timeout=timeout)
        elapsed = round(time.monotonic() - started, 3)
        log = output_dir / f"{name}.sh.log"
        log.write_text(completed.stdout +
                       (f"\n[stderr]\n{completed.stderr}" if completed.stderr else ""))
        checks[f"{name}.sh"] = {
            "returncode": completed.returncode,
            "script_sha256": sha256_file(script),
            "log_sha256": sha256_file(log),
            "runtime_seconds": elapsed,
            "command": ["bash", f"script/{name}.sh"],
        }
        if completed.returncode != 0:
            raise AcceptanceError(
                f"health check failed: {name}.sh exited {completed.returncode}")
    value = {"schema": HEALTH_SCHEMA,
             "commit": commit if commit is not None else git_revs(root)["head"],
             "checks": checks}
    write_json(output_dir / "health.json", value)
    return value


def read_health(root, health_dir, *, phase, inventory) -> dict:
    health_dir = Path(health_dir)
    path = health_dir / "health.json"
    if not path.is_file():
        if phase == "final":
            raise MissingArtifact(f"missing health evidence: {rel(root, path)}")
        return {"schema": HEALTH_SCHEMA, "checks": {}, "verify_all": VERIFY_ALL_NOTE,
                "status": "pending"}
    value = load_json(path)
    if value.get("schema") != HEALTH_SCHEMA:
        raise AcceptanceError(f"health schema mismatch: {rel(root, path)}")
    checks = {}
    for name in ("check_paths.sh", "check_tracked.sh", "verify_all.sh"):
        entry = (value.get("checks") or {}).get(name)
        if entry is None:
            continue
        log = health_dir / f"{name}.log"
        if not log.is_file():
            raise MissingArtifact(f"health transcript missing: {rel(root, log)}")
        if sha256_file(log) != entry.get("log_sha256"):
            raise DigestMismatch(f"health transcript digest mismatch: {rel(root, log)}")
        script = Path(root) / "script" / name
        if not script.is_file():
            raise MissingArtifact(f"health script missing: script/{name}")
        if sha256_file(script) != entry.get("script_sha256"):
            raise DigestMismatch(f"health script digest mismatch: script/{name}")
        if entry.get("returncode") != 0:
            raise AcceptanceError(f"health check failed: {name} exited "
                                  f"{entry.get('returncode')}")
        _pin_health(inventory, root, log)
        _pin(inventory, "health", root, script)
        checks[name] = {
            "returncode": entry.get("returncode"),
            "script_sha256": entry.get("script_sha256"),
            "log": _label_for(root, log), "log_sha256": entry.get("log_sha256"),
            "runtime_seconds": entry.get("runtime_seconds"),
            "command": entry.get("command", ["bash", f"script/{name}"]),
        }
    for name in ("check_paths.sh", "check_tracked.sh"):
        if name not in checks:
            raise MissingArtifact(f"missing health transcript for {name}")
    _pin_health(inventory, root, path)
    if "verify_all.sh" in checks:
        verify_all = {"status": "passed", **checks["verify_all.sh"]}
    elif phase == "final":
        raise AcceptanceError("verify_all evidence is required for a sealed receipt: run "
                              "`bash script/verify_all.sh` alone, capture its transcript in "
                              "health.json, then rebuild")
    else:
        verify_all = dict(VERIFY_ALL_NOTE)
    return {"schema": HEALTH_SCHEMA, "commit": value.get("commit"), "checks": checks,
            "verify_all": verify_all, "health_json_sha256": sha256_file(path)}


# --------------------------------------------------------------------------- #
# accepted inputs
# --------------------------------------------------------------------------- #
def verify_inputs(root, inputs_root, freeze_products, *, inventory):
    release = Path(inputs_root) / "manifest.json"
    if not release.is_file():
        raise MissingArtifact(f"missing accepted 5k release: {rel(root, release)}")
    release_digest = sha256_file(release)
    _pin(inventory, "inputs", root, release)
    axis = Path(root) / "stage2/01_bases/output/gene-axis/gene_axis_5000.csv"
    sealed = {}
    for source in SOURCES:
        sample = next(product for key, product in freeze_products.items()
                      if key.startswith(f"{source}-"))
        sealed[source] = {"dir": f"stage2/04_tsg_dose/execution/checkpoints/"
                                 f"round1-axis-local-inputs/{source}",
                          "digests": sample["sealed"]}
    return {
        "view_5k_release": {"path": rel(root, Path(inputs_root)),
                            "manifest_sha256": release_digest,
                            "manifest_bytes": release.stat().st_size},
        "gene_axis": {"path": rel(root, axis), "sha256": sha256_file(axis)}
        if axis.is_file() else {"path": rel(root, axis), "sha256": None},
        "sealed_step04_inputs": sealed,
    }


# --------------------------------------------------------------------------- #
# Steps 01-06 invariance proof
# --------------------------------------------------------------------------- #
def _git(root, *arguments, binary=False):
    return subprocess.run(["git", "-C", str(root), *arguments], capture_output=True,
                          text=not binary)


def _tracked_paths(root, step, rev=None):
    if rev is None:
        completed = _git(root, "ls-files", "-z", "--", f"stage2/{step}", binary=True)
        if completed.returncode == 0:
            return sorted(item.decode() for item in completed.stdout.split(b"\0") if item)
        base = Path(root) / "stage2" / step
        return sorted(path.relative_to(root).as_posix()
                      for path in base.rglob("*") if path.is_file())
    completed = _git(root, "ls-tree", "-r", "-z", rev, "--", f"stage2/{step}",
                     binary=True)
    if completed.returncode != 0:
        return None
    entries = []
    for raw in completed.stdout.split(b"\0"):
        if not raw:
            continue
        meta, _, path = raw.partition(b"\t")
        entries.append((meta.split()[2].decode(), path.decode()))
    return entries


def _worktree_digests(root, paths):
    records = {}
    for path in paths:
        full = Path(root) / path
        if full.is_file():
            records[path] = {"sha256": sha256_file(full), "bytes": full.stat().st_size}
    return records


def _rev_digests(root, step, rev):
    entries = _tracked_paths(root, step, rev=rev)
    if entries is None:
        return None
    completed = subprocess.run(
        ["git", "-C", str(root), "cat-file", "--batch"],
        input=("\n".join(blob for blob, _ in entries) + "\n").encode(),
        capture_output=True)
    if completed.returncode != 0:
        return None
    data = completed.stdout
    position, records = 0, {}
    for blob, path in entries:
        newline = data.index(b"\n", position)
        header = data[position:newline].split()
        if len(header) < 3 or header[1] != b"blob":
            return None
        size = int(header[2])
        body = data[newline + 1:newline + 1 + size]
        position = newline + 1 + size + 1
        records[path] = {"sha256": hashlib.sha256(body).hexdigest(), "bytes": size}
    return records


def _step_block(records):
    lines = "".join(f"{records[path]['sha256']}  {path}\n" for path in sorted(records))
    return sha256_text(lines), sum(records[path]["bytes"] for path in records)


def _sealed_artifacts(root, step):
    base = Path(root) / "stage2" / step
    paths = set()
    for pattern in ("output/**/manifest.json", "log/**/*acceptance*.json",
                    "log/**/health.json"):
        paths.update(path for path in base.glob(pattern) if path.is_file())
    return sorted(path for path in paths)


def _git_diff_paths(root, rev, step, *, branch="dev"):
    """Committed paths changed between the pinned planning ref and ``branch``."""
    for target in (branch, "HEAD"):
        completed = _git(root, "diff", "--name-only", f"{rev}..{target}", "--",
                         f"stage2/{step}")
        if completed.returncode == 0:
            return sorted(path for path in completed.stdout.splitlines() if path.strip()), target
    return None, None


def _worktree_paths(root, step):
    """Tracked paths with uncommitted changes (staged, unstaged or untracked)."""
    completed = _git(root, "status", "--porcelain", "--untracked-files=all", "--",
                     f"stage2/{step}")
    if completed.returncode != 0:
        return None
    return sorted(line[3:] for line in completed.stdout.splitlines() if line.strip())


def _git_diff_evidence(root, rev, diff_paths, dirty_paths, steps, branch="dev"):
    """Summarise the committed and working-tree evidence for the receipt."""
    committed = sorted({path for values in diff_paths.values()
                        for path in (values or [])})
    dirty = sorted({path for values in dirty_paths.values()
                    for path in (values or [])})
    unavailable = rev is None or any(value is None for value in diff_paths.values()) \
        or any(value is None for value in dirty_paths.values())
    if unavailable:
        status = "unavailable"
    else:
        status = "PASS" if not committed and not dirty else "CHANGED"
    command = None if rev is None else [
        "git", "diff", "--name-only", f"{rev}..{branch}", "--",
        *[f"stage2/{step}" for step in steps]]
    return {"rev": rev, "branch": branch, "command": command,
            "paths": committed, "empty": (None if unavailable else not committed),
            "status": status,
            "worktree_clean": (None if unavailable else not dirty),
            "worktree_paths": dirty,
            "per_step": {step: list(values or []) for step, values in diff_paths.items()},
            "per_step_worktree": {step: list(values or [])
                                  for step, values in dirty_paths.items()}}


def steps_invariance(root, baseline_path=None, refs_path=None, *, steps=STEPS,
                     prior_receipt=None):
    """Prove Steps 01-06 are unchanged.

    The reference is either a per-path manifest (``sha256sum`` lines) or the opaque
    single digest left by the planning capture. A per-path manifest pins the *on-disk*
    artifact inventory -- including gitignored payloads such as checkpoints and output
    manifests that ``git ls-files`` does not see -- so it is compared directly against
    the current files. The pinned planning commit and the dirty working tree are checked
    alongside it, so a source-code edit the artifact inventory does not list still fails.
    """
    root = Path(root)
    paths = default_paths(root)
    baseline_path = Path(baseline_path) if baseline_path else paths["baseline"]
    refs_path = Path(refs_path) if refs_path else paths["refs"]
    if not baseline_path.is_file():
        raise MissingArtifact(f"missing Steps 01-06 baseline: {rel(root, baseline_path)}")
    text = baseline_path.read_text()
    manifest, opaque = {}, None
    for line in text.splitlines():
        match = _BASELINE_LINE.match(line.strip())
        if match and match.group(2) != "-":
            manifest[match.group(2)] = match.group(1)
        elif match and match.group(2) == "-":
            opaque = match.group(1)
    mode = "manifest" if manifest else "opaque"
    rev = None
    if refs_path.is_file():
        tokens = refs_path.read_text().split()
        if tokens and re.fullmatch(r"[0-9a-f]{40}", tokens[0]):
            rev = tokens[0]
    current, current_blocks = {}, {}
    for step in steps:
        if manifest:
            # The manifest is the declared artifact inventory: verify exactly those
            # paths on disk, not the (smaller) git-tracked subset.
            pinned = sorted(path for path in manifest
                            if path.startswith(f"stage2/{step}/"))
            current[step] = _worktree_digests(root, pinned)
        else:
            current[step] = _worktree_digests(root, _tracked_paths(root, step))
        current_blocks[step] = _step_block(current[step])
    reference = None
    if rev is not None and not manifest:
        reference = {}
        for step in steps:
            records = _rev_digests(root, step, rev)
            if records is None:
                reference = None
                break
            reference[step] = records
    combined = sha256_text("".join(f"{current_blocks[step][0]}  stage2/{step}\n"
                                   for step in steps))
    matches_recipe = bool(opaque and combined == opaque)
    if not manifest and reference is None:
        raise InvarianceChanged(
            "invariance baseline is opaque and the pinned planning ref is unavailable")
    if not manifest and not matches_recipe and reference is None:
        raise InvarianceChanged("invariance baseline cannot be reproduced")
    diff_paths, diff_branch = {}, None
    if rev is not None:
        for step in steps:
            changed, target = _git_diff_paths(root, rev, step)
            diff_paths[step] = changed
            diff_branch = diff_branch or target
    dirty_paths = {step: _worktree_paths(root, step) for step in steps}
    proof_steps, verdict = {}, "PASS"
    for step in steps:
        if manifest:
            baseline_records = {path: digest for path, digest in manifest.items()
                                if path.startswith(f"stage2/{step}/")}
        else:
            baseline_records = reference[step]
        normalized = {}
        for path, value in (baseline_records or {}).items():
            if isinstance(value, dict):
                normalized[path] = value
            else:
                normalized[path] = {"sha256": value,
                                    "bytes": current[step].get(path, {}).get("bytes", 0)}
        baseline_digest, baseline_bytes = (_step_block(normalized) if normalized
                                           else (sha256_text(""), 0))
        changed = []
        if manifest:
            # Additions cannot hide: current is read from the pinned path set, so a
            # missing path is a change and every present path is compared one-to-one.
            for path in sorted(normalized):
                before = normalized[path]["sha256"]
                after = current[step].get(path, {}).get("sha256")
                if before != after:
                    changed.append({"path": path, "baseline": before, "current": after})
        else:
            for path in sorted(set(normalized) | set(current[step])):
                before = normalized.get(path, {}).get("sha256")
                after = current[step].get(path, {}).get("sha256")
                if before != after:
                    changed.append({"path": path, "baseline": before, "current": after})
        committed = diff_paths.get(step)
        dirty = dirty_paths.get(step)
        unavailable = rev is not None and (committed is None or dirty is None)
        step_verdict = ("PASS" if not changed and not committed and not dirty
                        and not unavailable else "CHANGED")
        if step_verdict != "PASS":
            verdict = "CHANGED"
        proof_steps[step] = {
            "verdict": step_verdict,
            "files": len(current[step]),
            "bytes": current_blocks[step][1],
            "digest": current_blocks[step][0],
            "baseline_digest": baseline_digest,
            "baseline_files": len(baseline_records),
            "baseline_bytes": baseline_bytes,
            "changed": changed,
            "committed_changes": committed,
            "worktree_changes": dirty,
        }
    sealed = {}
    for step in steps:
        records = _worktree_digests(root, [rel(root, path)
                                           for path in _sealed_artifacts(root, step)])
        digest, total = _step_block(records) if records else (sha256_text(""), 0)
        sealed[step] = {"files": len(records), "digest": digest, "bytes": total}
    sealed_blob = sha256_text("".join(
        f"{sealed[step]['digest']}  stage2/{step}\n" for step in steps))
    if prior_receipt is not None:
        prior = load_json(prior_receipt)
        block = (prior.get("steps01_06") or {})
        for step in steps:
            before = (block.get("steps") or {}).get(step, {}).get("digest")
            if before is not None and before != proof_steps[step]["digest"]:
                verdict = "CHANGED"
                proof_steps[step]["verdict"] = "CHANGED"
                proof_steps[step]["changed"] = [{"path": f"stage2/{step}", "baseline":
                                                 before, "current":
                                                 proof_steps[step]["digest"]}]
        if block.get("sealed", {}).get("digest") not in (None, sealed_blob):
            verdict = "CHANGED"
    value = {
        "schema": INVARIANCE_SCHEMA,
        "verdict": verdict,
        "steps": proof_steps,
        "sealed": {"files": sum(item["files"] for item in sealed.values()),
                   "digest": sealed_blob, "per_step": sealed},
        "combined_digest": combined,
        "baseline_file": {"path": rel(root, baseline_path),
                          "sha256": sha256_file(baseline_path), "mode": mode,
                          "matches_recipe": matches_recipe,
                          "opaque_digest": opaque,
                          "files": len(manifest),
                          "lines": len(text.splitlines())},
        "baseline_ref": {"rev": rev, "source": rel(root, refs_path) if refs_path.is_file()
                         else None,
                         "refs_sha256": sha256_file(refs_path)
                         if refs_path.is_file() else None},
        "git_diff": _git_diff_evidence(root, rev, diff_paths, dirty_paths, steps,
                                       diff_branch or "dev"),
    }
    if prior_receipt is not None and verdict != "PASS":
        raise InvarianceChanged("Steps 01-06 changed since the recorded baseline receipt")
    return value


# --------------------------------------------------------------------------- #
# receipt assembly
# --------------------------------------------------------------------------- #
def _portable_argv(argv, root):
    root = Path(root).resolve()
    portable = []
    for argument in argv:
        text = str(argument)
        candidate = Path(text)
        if candidate.is_absolute():
            try:
                portable.append("./" + candidate.resolve().relative_to(root).as_posix())
                continue
            except ValueError:
                portable.append(candidate.name)
                continue
        portable.append(text)
    return portable


def build_receipt(root, *, phase="partial", timestamp=None, environment=None,
                  commands=None, views=None, freezes=None, truth=None, runs=None,
                  report=None, figure=None, log=None, inputs=None, health_dir=None,
                  baseline=None, refs=None, baseline_receipt=None, preregistration=None,
                  preregistration_text=None, verify_binary=False, python=None,
                  replay_figures=True):
    if phase not in ("partial", "final"):
        raise ValueError("phase must be 'partial' or 'final'")
    root = Path(root)
    paths = default_paths(root)
    views_root = Path(views) if views else paths["views"]
    freezes_root = Path(freezes) if freezes else paths["freezes"]
    truth_root = Path(truth) if truth else paths["truth"]
    runs_root = Path(runs) if runs else paths["runs"]
    report_root = Path(report) if report else paths["report"]
    figure_root = Path(figure) if figure else paths["figure"]
    log_root = Path(log) if log else paths["log"]
    health_root = Path(health_dir) if health_dir else log_root / "acceptance/health"
    inputs_root = Path(inputs) if inputs else paths["inputs"]
    baseline_path = Path(baseline) if baseline else paths["baseline"]
    refs_path = Path(refs) if refs else paths["refs"]
    if phase == "final" and not verify_binary:
        raise AcceptanceError("a sealed receipt requires --verify-binary")
    inventory = {}
    view_products = verify_views(root, views_root, phase=phase,
                                 verify_binary=verify_binary, inventory=inventory)
    freeze_products = verify_freezes(root, freezes_root, view_products, phase=phase,
                                     verify_binary=verify_binary, inventory=inventory)
    truth_products = verify_truth(root, truth_root, phase=phase,
                                  verify_binary=verify_binary, inventory=inventory)
    runs_value = verify_runs(root, runs_root, view_products, freeze_products,
                             truth_products, phase=phase,
                             verify_binary=verify_binary, inventory=inventory)
    smoke = verify_smoke(root, paths["step"], phase=phase, inventory=inventory)
    report_value = verify_report(root, report_root, figure_root, phase=phase,
                                 inventory=inventory, python=python,
                                 replay_figures=replay_figures)
    health = read_health(root, health_root, phase=phase, inventory=inventory)
    proof = steps_invariance(root, baseline_path, refs_path,
                             prior_receipt=baseline_receipt)
    inputs_value = verify_inputs(root, inputs_root, freeze_products, inventory=inventory)
    if preregistration_text is not None:
        preregistration_sha256 = sha256_text(preregistration_text)
    elif preregistration is not None:
        preregistration_sha256 = sha256_file(preregistration)
    else:
        if phase == "final":
            raise AcceptanceError("pre-registration digest required for a sealed receipt")
        preregistration_sha256 = None
    environment = environment if environment is not None else environment_probe(python)
    if commands is None:
        commands = []
        for name, entry in health.get("checks", {}).items():
            commands.append({"name": name, "argv": entry.get("command", [])})
    if phase == "final" and proof["verdict"] != "PASS":
        raise InvarianceChanged("Steps 01-06 changed: the sealed receipt cannot be built")
    pending = runs_value["pending_origins"]
    report_scripts = {}
    for name in ("execution/report.py", "execution/render.py", "figure/script/render.py"):
        path = paths["step"] / name
        report_scripts[name] = sha256_file(path) if path.is_file() else None
        if phase == "final" and report_scripts[name] is None:
            raise MissingArtifact(f"missing report script: stage2/tsg_view/{name}")
    replay = {"smoke": smoke, "matrix": [runs_value["replay"][origin_key(entry["origin"])]
                                         for entry in runs_value["present"]],
              "all_exact": bool(smoke) and all(
                  entry.get("validation_replay") == "exact" and
                  entry.get("test_replay") == "exact" for entry in
                  runs_value["replay"].values())}
    receipt = {
        "schema": SCHEMA,
        "status": "passed" if phase == "final" else "partial",
        "phase": phase,
        "generated_at": timestamp or now_utc(),
        "study": STUDY,
        "scenario": SCENARIO,
        "seeds": {"main": list(SEEDS), "solo_deepest": list(SOLO_DEEPEST_SEEDS),
                  "solo_draw": SOLO_SEED},
        "formal_fits": EXPECTED_FITS,
        "complete_origins": len(runs_value["present"]),
        "pending_origins": runs_value["pending_origins"],
        "plan_sha256": runs_value["plan_sha256"],
        "plan_manifest_sha256": runs_value["plan_manifest_sha256"],
        "input_manifest_sha256": inputs_value["view_5k_release"]["manifest_sha256"],
        "inputs": inputs_value,
        "environment": environment,
        "caps": {"memory_max": CAP_TEXT[0], "memory_swap_max": CAP_TEXT[1],
                 "memory_max_bytes": MEMORY_MAX, "swap_max_bytes": SWAP_MAX,
                 "telemetry_seconds": TELEMETRY_SECONDS},
        "inventory": {group: dict(sorted(records.items()))
                      for group, records in sorted(inventory.items())},
        "products": {"views": view_products, "freezes": freeze_products,
                     "truth": truth_products, "report": report_value,
                     "figures": report_value["figures"]},
        "origins": runs_value["present"],
        "resources": dict(sorted(runs_value["resources"].items())),
        "replay": replay,
        "formal_cap": runs_value["formal_cap"],
        "run_records": runs_value["run_records"],
        "formal_verify": runs_value["formal_verify"],
        "report": report_value,
        "figures": {"render_status_sha256":
                    report_value["csv_only_replay"].get("render_status_sha256"),
                    "files": report_value["figures"],
                    "suite": report_value["csv_only_replay"].get("suite")},
        "csv_only_replay": report_value["csv_only_replay"],
        "report_scripts_sha256": report_scripts,
        "steps01_06": proof,
        "health": health,
        "commands": commands,
        "git": git_revs(root),
        "pre_registration_sha256": preregistration_sha256,
        "pending": [] if phase == "final" else pending,
    }
    if phase == "final":
        for key in ("plan_sha256", "input_manifest_sha256"):
            if not receipt[key]:
                raise MissingArtifact(f"missing {key}")
    return receipt


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subs = parser.add_subparsers(dest="action", required=True)
    health = subs.add_parser("health", help="run the fast offline health checks")
    health.add_argument("--root", type=Path, default=Path("."))
    health.add_argument("--output", type=Path, required=True)
    health.add_argument("--verify-all", action="store_true",
                        help="also run the heavy script/verify_all.sh, alone")
    health.add_argument("--timeout", type=float, default=900.0)
    receipt = subs.add_parser("receipt", help="assemble the acceptance receipt")
    receipt.add_argument("--root", type=Path, default=Path("."))
    receipt.add_argument("--phase", choices=("partial", "final"), default="partial")
    receipt.add_argument("--views", type=Path)
    receipt.add_argument("--freezes", type=Path)
    receipt.add_argument("--truth", type=Path)
    receipt.add_argument("--runs", type=Path)
    receipt.add_argument("--report", type=Path)
    receipt.add_argument("--figure", type=Path)
    receipt.add_argument("--log", type=Path)
    receipt.add_argument("--inputs", type=Path)
    receipt.add_argument("--health-dir", type=Path)
    receipt.add_argument("--output", type=Path)
    receipt.add_argument("--baseline", type=Path)
    receipt.add_argument("--refs", type=Path)
    receipt.add_argument("--baseline-receipt", type=Path)
    receipt.add_argument("--preregistration", type=Path)
    receipt.add_argument("--verify-binary", action="store_true")
    receipt.add_argument("--no-figure-replay", action="store_true")
    receipt.add_argument("--timestamp")
    return parser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = _parser()
    arguments = parser.parse_args(argv)
    if arguments.action == "health":
        value = run_health(arguments.root, arguments.output,
                           verify_all=arguments.verify_all,
                           timeout=arguments.timeout)
        sys.stdout.write(json.dumps(
            {"schema": value["schema"],
             "checks": {name: {"returncode": entry["returncode"],
                               "runtime_seconds": entry["runtime_seconds"]}
                        for name, entry in value["checks"].items()}},
            indent=2, sort_keys=True) + "\n")
        return 0
    paths = default_paths(arguments.root)
    output = arguments.output or paths["output"]
    receipt = build_receipt(
        arguments.root, phase=arguments.phase,
        views=arguments.views, freezes=arguments.freezes, truth=arguments.truth,
        runs=arguments.runs, report=arguments.report, figure=arguments.figure,
        log=arguments.log, inputs=arguments.inputs, health_dir=arguments.health_dir,
        baseline=arguments.baseline, refs=arguments.refs,
        baseline_receipt=arguments.baseline_receipt,
        preregistration=arguments.preregistration,
        verify_binary=arguments.verify_binary, timestamp=arguments.timestamp,
        replay_figures=not arguments.no_figure_replay,
        commands=[{"name": "receipt", "argv": _portable_argv(
            ["python", "stage2/tsg_view/execution/acceptance.py", *argv],
            arguments.root)}])
    write_receipt(receipt, output)
    sys.stdout.write(json.dumps({
        "receipt": str(output), "status": receipt["status"],
        "complete_origins": receipt["complete_origins"],
        "pending_origins": len(receipt["pending_origins"]),
        "steps01_06": receipt["steps01_06"]["verdict"],
        "csv_only_replay": receipt["csv_only_replay"].get("byte_identical"),
    }, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
