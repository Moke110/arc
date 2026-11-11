"""Re-hash the accepted bases payloads once, serially, under the measured caps (#156 E1).

This proves a path move did not rewrite bytes. It reads the sealed manifests at
their stable paths, checks each manifest's own digest against the accepted seal,
then hashes every file the manifest declares -- including the multi-gigabyte
``*.h5ad`` and ``*.npy`` payloads -- and compares it with the declared digest.
The rejected historical Tahoe training payload is deliberately skipped: it is
not an accepted dependency and #155 records it as rejected and superseded.

    systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
        env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
        /usr/bin/time -v python3 stage2/01_bases/execution/verify_payload_digests.py \
            --step-root /home/chang/arc/stage2/01_bases \
            --record /path/to/log/156-payload-digests.json
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

SCHEMA = "bases-156-payload-digests-v1"

# Accepted seal digests (155-migration-map.json acceptance gates and the #155
# metadata-binding audit). A moved manifest keeps its bytes, so these stay valid.
SEALS = {
    "output/gene-identity/manifest.json": "5f7517c34d7fda2ff61993b625086bd3a5db1de51c097b6cf9f87aec7bf98ade",
    "output/tenv-raw-counts/tenv_manifest.json": "e0941b240985ed0c77687d363c8e2ba1b069d86f1227b1d512afad6c4897aa5b",
    "output/tahoe-raw-counts/manifest.json": "7ea21f45e0a9d79383f243a5dd518a7e2851433da400e06e31e575a53cdada44",
    "output/shared-expression/manifest.json": "fb352f2be2c93408b58daddf5f0ed9d2972b97c0f9a5f16b187316c6245987c3",
    "output/training-corpora/manifest.json": "aac36d2d8306c234da581ae2763b8703a323db89870c7db42a15fe0fd5be564b",
    "output/tenv-truth/manifest.json": "ab73b3018bfaaaedb91328f455c1f1bf3832949de553026923c06d4d11f8977f",
    "output/tenv-hvg/manifest.json": "953103a1093e4364158075553dcb036c560c14fa61ea95fa8f4157a8a881c3a6",
    "output/historical-training-release/manifest.json": "3543df76e728b3d25d63cda20e3187e96fbca96ac370d38e634e25405835a1a4",
}
# The accepted plotted values and figure before the extension (audit evidence copies).
ACCEPTED_EVIDENCE = {
    "log/diagnostic/diagnostic_pr.csv": "eaad870b346ead01f08bb4af17565e8f62a250168dbc70c4e0d984b52281ac1a",
    "log/diagnostic/diagnostic_pr.png": "35d14266283b296c1585af3d58f0bacb15a3d7f854cde4972c433375bbdd911f",
}
SKIP_DECLARED = {"output/historical-training-release/tahoe_bulk.h5ad"}


def digest(path):
    value = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            value.update(chunk)
            size += len(chunk)
    return value.hexdigest(), size


def resolve(step, product, name):
    if name == "hgnc_complete_set.txt":
        return step.parents[1] / "pert_data/Annotation/source/hgnc_complete_set.txt"
    if name in ("tahoe_mapping.csv", "tenv_mapping.csv"):
        return step / "log/historical/shared-v1-mappings" / name
    if product.endswith("gene-identity") and name in ("genes.csv", "gene_crosswalk.csv"):
        return step / "output/gene-identity/tables" / name
    return step / product / name


def declared(step, product):
    """Return (relative_path, expected_digest) pairs the product manifest declares."""
    root = step / product
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        manifest_path = root / "tenv_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    pairs = []  # (absolute path, declared digest)
    if "genes_sha256" in manifest:
        pairs.append((root / "genes.csv", manifest["genes_sha256"]))
    if "files" in manifest and isinstance(manifest["files"], dict):
        for name, value in manifest["files"].items():
            pairs.append((root / name, value))
    if "mapping_sha256" in manifest:
        pairs.append((root / "token_mapping.csv", manifest["mapping_sha256"]))
    if "h5ad_sha256" in manifest:
        pairs.append((root / "tahoe.h5ad", manifest["h5ad_sha256"]))
    if "cohorts" in manifest:
        pairs.append((root / manifest["cohorts"]["file"], manifest["cohorts"]["sha256"]))
    for population in manifest.get("populations", {}).values():
        if isinstance(population, dict) and "ranking" in population:
            pairs.append((root / population["ranking"], population["sha256"]))
    if "inputs" in manifest and isinstance(manifest["inputs"], dict) and "genes_sha256" in manifest:
        for name, value in manifest["inputs"].items():
            pairs.append((resolve(step, product, name), value["sha256"]))
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step-root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--record", type=Path)
    parser.add_argument("--fast", action="store_true",
                        help="skip .h5ad/.npy files for a metadata-only re-check")
    args = parser.parse_args()
    step = args.step_root.resolve()
    started = time.time()
    checks = []
    failures = []
    total_bytes = 0
    for relative, expected in SEALS.items():
        path = step / relative
        if not path.is_file():
            checks.append({"path": relative, "status": "MISSING"})
            failures.append(relative)
            continue
        actual, size = digest(path)
        total_bytes += size
        status = "PASS" if actual == expected else "FAIL"
        checks.append({"path": relative, "declared": expected, "actual": actual,
                       "bytes": size, "role": "manifest seal", "status": status})
        if status != "PASS":
            failures.append(relative)
    for relative, prefix in ACCEPTED_EVIDENCE.items():
        path = step / relative
        if not path.is_file():
            checks.append({"path": relative, "status": "MISSING"})
            failures.append(relative)
            continue
        actual, size = digest(path)
        total_bytes += size
        status = "PASS" if actual.startswith(prefix) else "FAIL"
        checks.append({"path": relative, "declared_prefix": prefix, "actual": actual,
                       "bytes": size, "role": "accepted evidence", "status": status})
        if status != "PASS":
            failures.append(relative)
    payloads = set()
    for product in [item.rsplit("/", 1)[0] for item in SEALS]:
        for path, expected in declared(step, product):
            if str(path) in payloads:
                continue
            payloads.add(str(path))
            try:
                relative = path.relative_to(step)
            except ValueError:
                relative = path.relative_to(step.parents[1])
            if str(relative) in SKIP_DECLARED:
                checks.append({"path": str(relative), "status": "SKIPPED",
                               "reason": "rejected superseded historical Tahoe payload; not an input"})
                continue
            if args.fast and path.suffix in {".h5ad", ".npy", ".npz", ".bin"}:
                checks.append({"path": str(relative), "status": "SKIPPED_FAST"})
                continue
            if not path.is_file():
                checks.append({"path": str(relative), "status": "MISSING",
                               "declared": expected})
                failures.append(str(relative))
                continue
            actual, size = digest(path)
            total_bytes += size
            status = "PASS" if actual == expected else "FAIL"
            checks.append({"path": str(relative), "declared": expected, "actual": actual,
                           "bytes": size, "role": "declared payload", "status": status})
            if status != "PASS":
                failures.append(str(relative))
    successor = step / "output/tenv-truth-successor"
    if (successor / "manifest.json").is_file():
        manifest = json.loads((successor / "manifest.json").read_text())
        worker = successor / "deseq2_contrasts.R"
        actual, size = digest(worker)
        total_bytes += size
        status = "PASS" if actual == manifest["worker_sha256"] else "FAIL"
        checks.append({"path": "output/tenv-truth-successor/deseq2_contrasts.R",
                       "declared": manifest["worker_sha256"], "actual": actual,
                       "bytes": size, "role": "successor worker seal", "status": status})
        if status != "PASS":
            failures.append("output/tenv-truth-successor/deseq2_contrasts.R")
        checks.append({"path": "output/tenv-truth-successor/manifest.json",
                       "actual": digest(successor / "manifest.json")[0],
                       "role": "re-issued successor manifest (digest recorded by #156)",
                       "status": "RECORDED"})
    record = {"schema": SCHEMA, "ticket": 156, "parent": 154,
              "step_root": str(step), "expression_payloads_opened": True,
              "fast": args.fast, "command": [sys.executable, *sys.argv],
              "result": "PASS" if not failures else "FAIL",
              "failures": failures,
              "checked_files": len(checks), "bytes_hashed": total_bytes,
              "checks": checks,
              "seconds": round(time.time() - started, 2),
              "peak_rss_kib": int(open("/proc/self/status").read().split("VmHWM:")[1].split()[0])}
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: record[key] for key in
                      ("result", "checked_files", "bytes_hashed", "seconds", "peak_rss_kib")},
                     indent=2))
    if failures:
        print("FAILURES:", json.dumps(failures, indent=2))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
