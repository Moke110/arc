"""Re-issue the truth successor envelope at its stable path, never editing it (#156 E2).

The accepted envelope records absolute input paths that
``bases.truth_successor.check_successor()`` resolves and hashes at runtime, so a
path move requires a new envelope. This script:

1. refuses to touch a stable envelope that already verifies against the stable
   training release (idempotent re-run);
2. moves the superseded envelope, bytes intact, to
   ``log/historical/seals/truth-successor-envelope-v1/`` (nothing is deleted);
3. creates a candidate envelope from the new absolute paths,
4. fully re-verifies the candidate (the complete equivalence proof, including
   the original raw-count replay) -- a failed candidate leaves the accepted
   products untouched and is not published;
5. only after PASS swaps the candidate into ``output/tenv-truth-successor``.

Run it as one serial capped job::

    systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
        env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
        /usr/bin/time -v python3 stage2/01_bases/execution/reissue_truth_successor.py \
            --step-root /home/chang/arc/stage2/01_bases \
            --candidate /home/chang/arc/.scratch/stage2-bases/runtime/v3/successor-candidate \
            --record /path/to/log/156-successor-reissue.json
"""
import argparse
import hashlib
import json
from pathlib import Path
import os
import shutil
import subprocess
import sys
import time

import public_api  # noqa: F401  (binds the single public implementation in output/api)

SCHEMA = "bases-156-truth-successor-reissue-v1"
STEP_ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verified(envelope, training):
    """True when the envelope already binds this training release and verifies."""
    from bases.truth_successor import check_successor

    envelope, training = Path(envelope), Path(training)
    if not (envelope / "manifest.json").is_file():
        return False, None
    try:
        value, _ = check_successor(envelope, training, full_payload=False)
    except (ValueError, OSError, KeyError, FileNotFoundError):
        return False, None
    return True, value


def run(command, record, label, output_path=None):
    started = time.time()
    result = subprocess.run(command, capture_output=True, text=True)
    record.setdefault("commands", []).append(
        {"label": label, "argv": [str(part) for part in command],
         "returncode": result.returncode, "seconds": round(time.time() - started, 2),
         "stdout_tail": result.stdout.strip().splitlines()[-3:],
         "stderr_tail": result.stderr.strip().splitlines()[-3:]})
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed: {result.stderr.strip()[-800:]}")
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(result.stdout)
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step-root", type=Path, default=STEP_ROOT)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--record", type=Path)
    parser.add_argument("--records-dir", type=Path,
                        help="where the create/verify JSON results are recorded "
                             "(defaults to step-root/log/records)")
    parser.add_argument("--rscript", type=Path,
                        default=Path("/home/chang/miniconda3/envs/biomni_e1/bin/Rscript"))
    args = parser.parse_args()
    step = args.step_root.resolve()
    execution = step / "execution"
    output = step / "output"
    records = (args.records_dir.resolve() if args.records_dir else step / "log/records")
    seals = step / "log/historical/seals"
    stable = output / "tenv-truth-successor"
    training = output / "training-corpora"
    candidate = args.candidate.resolve()
    record = {"schema": SCHEMA, "ticket": 156, "parent": 154, "step_root": str(step),
              "command": [sys.executable, *sys.argv], "started": time.time()}
    if not training.is_dir() or not (records / "truth_complete_verify.json").is_file():
        print("ERROR: the stable training release or accepted verification record is missing",
              file=sys.stderr)
        return 2
    ok, value = verified(stable, training)
    if ok:
        record.update({"result": "ALREADY_REISSUED",
                       "manifest_sha256": digest(stable / "manifest.json"),
                       "bindings": {name: item["path"] for name, item in value["bindings"].items()}})
        print(json.dumps(record, indent=2, sort_keys=True))
        return 0
    if candidate.exists():
        shutil.rmtree(candidate)
    create = [sys.executable, str(execution / "bind_truth_successor.py"),
              "--truth", str(output / "tenv-truth"),
              "--original-training", str(output / "historical-training-release"),
              "--training", str(training),
              "--original-native", str(output / "tenv-raw-counts/tenv.h5ad"),
              "--native", str(output / "tenv-raw-counts/tenv.h5ad"),
              "--accepted-verification", str(records / "truth_complete_verify.json"),
              "--output", str(candidate), "--rscript", str(args.rscript)]
    verify = [sys.executable, str(execution / "bind_truth_successor.py"), "--verify",
              "--training", str(training), "--output", str(candidate),
              "--rscript", str(args.rscript)]
    try:
        run(create, record, "create candidate",
            records / "156_truth_successor_build.json")
        candidate_manifest = digest(candidate / "manifest.json")
        run(verify, record, "verify candidate",
            records / "156_truth_successor_verify.json")
        record["candidate"] = {"path": str(candidate), "manifest_sha256": candidate_manifest,
                               "create_pass": True, "verify_pass": True}
    except Exception as error:  # noqa: BLE001 - a failed candidate must not publish
        record.update({"result": "CANDIDATE_FAILED", "problem": str(error),
                       "stable_untouched": digest(stable / "manifest.json")
                       if (stable / "manifest.json").is_file() else None})
        if args.record:
            args.record.parent.mkdir(parents=True, exist_ok=True)
            args.record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"result": "CANDIDATE_FAILED", "problem": str(error)}, indent=2))
        return 1
    old = digest(stable / "manifest.json") if (stable / "manifest.json").is_file() else None
    seals.mkdir(parents=True, exist_ok=True)
    preserved = seals / "truth-successor-envelope-v1"
    if preserved.exists():
        if (digest(preserved / "manifest.json") != old
                or digest(preserved / "deseq2_contrasts.R")
                != digest(stable / "deseq2_contrasts.R")):
            print("ERROR: a different superseded envelope is already preserved", file=sys.stderr)
            return 1
        shutil.rmtree(stable)
    else:
        if old is None:
            print("ERROR: no superseded envelope to preserve", file=sys.stderr)
            return 1
        os.rename(stable, preserved)
    os.rename(candidate, stable)
    record.update({"result": "PASS", "superseded_manifest_sha256": old,
                   "superseded_bytes_preserved_at": str(preserved),
                   "manifest_sha256": digest(stable / "manifest.json"),
                   "bindings": {name: item["path"] for name, item in
                                json.loads((stable / "manifest.json").read_text())["bindings"].items()}})
    record["finished"] = time.time()
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: record[key] for key in
                      ("result", "superseded_manifest_sha256", "manifest_sha256")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
