"""Bounded access to the frozen corpus, fitted artifacts and reconstruction scores."""

from contextlib import contextmanager
import csv
import hashlib
import json
from math import isclose
from numbers import Integral
import os
from pathlib import Path
import tempfile

import numpy as np


ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
STEP = ROOT / "stage2" / "07_autoencoder"
CORPUS_STEP = ROOT / "stage2" / "02_dataset_corpus"
DEFAULT_RELEASE = CORPUS_STEP / "output" / "corpus" / "shared-cm-v1"
DEFAULT_SPLIT = STEP / "output" / "split.json"
DEFAULT_MEAN = STEP / "output" / "models" / "mean"
AAE_MODEL_ID = "AAE_512"
AAE_TRANSFORMER_MODEL_ID = "AAE_512_Transformer"
AAE_NUISANCE_MODEL_ID = "AAE_512+16"
AAE_TRANSFORMER_MODEL_IDS = (AAE_TRANSFORMER_MODEL_ID, AAE_NUISANCE_MODEL_ID)
AAE_MODEL_IDS = (AAE_MODEL_ID, *AAE_TRANSFORMER_MODEL_IDS)
AAE_RECONSTRUCTION_TOKENS = 512
AAE_NUISANCE_TOKENS = 16
AAE_WIDTH = 128
MATRIX_AE_MODEL_ID = "MatrixAE_64x8"
MATRIX_AE_TOKENS = 64
MATRIX_AE_CHANNELS = 8
# A batch at this limit occupies 71 MiB of float32 expression on the real gene axis.
MAX_BATCH_SIZE = 1024
SOURCES = ("tahoe", "tenv")
ROLES = ("train", "validation", "test")
SEED = 20261109
EXPECTED_COUNTS = {
    "tahoe": ((10, 47822), (2, 9600), (2, 9596)),
    "tenv": ((10, 2810), (3, 843), (3, 825)),
}
REPLAY_TOLERANCE = {"rtol": 1e-6, "atol": 1e-8}


def aae_encoder_tokens(model_id):
    if model_id == AAE_NUISANCE_MODEL_ID:
        return AAE_RECONSTRUCTION_TOKENS + AAE_NUISANCE_TOKENS
    return AAE_RECONSTRUCTION_TOKENS


def file_sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def check_checkpoint_selection(history, metadata, model_id):
    json_bytes(history)
    best = min(history, key=lambda check: (check["validation_mse"]["tenv"], check["epoch"]))
    if (metadata["selected_epoch"] != best["epoch"]
            or metadata["selected_validation_mse"] != best["validation_mse"]):
        raise ValueError(f"{model_id}: checkpoint differs from earliest tenv minimum")


def validation_mse_matches(sources, selected_mse):
    return all(
        isclose(sources[source]["mse"], selected_mse[source],
                rel_tol=REPLAY_TOLERANCE["rtol"], abs_tol=REPLAY_TOLERANCE["atol"])
        for source in SOURCES
    )


def digest(value):
    return hashlib.sha256(json_bytes(value)).hexdigest()


def publish_json(path, value):
    """Publish a complete JSON record once, with an exclusive destination."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent) as staging:
        pending = Path(staging) / path.name
        pending.write_bytes(json_bytes(value))
        os.link(pending, path)


def csv_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


class Corpus:
    def __init__(self, release, fixture=False):
        release = Path(release)
        manifest_path = release / "manifest.json"
        manifest_bytes = manifest_path.read_bytes()
        self.manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
        manifest = json.loads(manifest_bytes)
        self.fixture = fixture
        self.normalization = manifest["normalization"]
        if (manifest["schema_version"] != 1 or manifest["fixture"] is not fixture
                or self.normalization != "log1p(CP10K of raw counts projected to the ordered CM axis)"):
            raise ValueError("unexpected Step 02 release schema, fixture mode or expression scale")
        if not fixture:
            verification_path = CORPUS_STEP / "execution" / "build_record" / "shared_cm_v1.json"
            verification = json.loads(verification_path.read_text())["verify"]
            if (verification["result"] != "PASS"
                    or verification["manifest_sha256"] != self.manifest_sha256):
                raise ValueError("corpus manifest differs from the verified Step 02 release")
        payloads = {"genes.csv", *(f"{source}_{suffix}" for source in SOURCES
                                  for suffix in ("expression.npy", "profiles.csv", "mapping.csv"))}
        files = manifest["files"]
        if set(files) != payloads or files["genes.csv"] != manifest["inputs"]["genes.csv"]:
            raise ValueError("unexpected Step 02 payload list or gene order")
        for name, sha256 in files.items():
            if file_sha256(release / name) != sha256:
                raise ValueError(f"corpus file digest mismatch: {name}")

        axis = [row["hgnc_id"] for row in csv_rows(release / "genes.csv")]
        if not axis or "" in axis or len(set(axis)) != len(axis):
            raise ValueError("invalid shared gene axis")
        self.gene_count = len(axis)
        self.gene_order_sha256 = files["genes.csv"]
        self.profiles = {}
        self._expression = {}
        for source in SOURCES:
            profiles = csv_rows(release / f"{source}_profiles.csv")
            identities = [row["profile_id"] for row in profiles]
            source_rows = [int(row["source_row"]) for row in profiles]
            if ("" in identities or len(set(identities)) != len(profiles)
                    or source_rows != sorted(set(source_rows)) or any(row < 0 for row in source_rows)
                    or any(not row["plate"] or row["is_control"] not in {"True", "False"}
                           for row in profiles)):
                raise ValueError(f"{source}: invalid profile identities or metadata")
            expression = np.load(release / f"{source}_expression.npy", mmap_mode="r")
            shape = (len(profiles), self.gene_count)
            if (expression.dtype != np.float32 or expression.shape != shape
                    or manifest["sources"][source]["shape"] != list(shape)):
                raise ValueError(f"{source}: unexpected expression shape or dtype")
            self.profiles[source] = profiles
            self._expression[source] = expression

    def read(self, source, rows):
        expression = self._expression[source]
        if isinstance(rows, Integral) and not isinstance(rows, bool):
            return expression[rows]
        sliced = isinstance(rows, slice)
        count = len(range(*rows.indices(len(expression)))) if sliced else len(rows)
        if count > MAX_BATCH_SIZE:
            raise ValueError(f"corpus reads are limited to at most {MAX_BATCH_SIZE} rows")
        if sliced:
            return expression[rows]
        indices = np.asarray(rows)
        if indices.ndim != 1 or (count and indices.dtype.kind not in "biu"):
            raise ValueError("corpus reads require one-dimensional row indices")
        return expression[indices if count else []]


def plate_rank(source, plate):
    # Preserve the original Step 07 serialization, including json.dumps spacing.
    return hashlib.sha256(json.dumps([SEED, source, plate]).encode()).digest()


def split_memberships(corpus):
    sources = {}
    rows = {}
    for source in SOURCES:
        profiles = corpus.profiles[source]
        plates = sorted({row["plate"] for row in profiles}, key=lambda plate: plate_rank(source, plate))
        held_out = max(1, len(plates) // 5)
        if len(plates) <= 2 * held_out:
            raise ValueError(f"{source}: split requires nonempty train, validation and test")
        groups = dict(train=plates[2 * held_out:], validation=plates[:held_out],
                      test=plates[held_out:2 * held_out])
        rows[source] = {}
        roles = {}
        for role, role_plates in groups.items():
            plate_set = set(role_plates)
            indices = np.array([i for i, row in enumerate(profiles) if row["plate"] in plate_set],
                               dtype=np.int64)
            indices.setflags(write=False)
            rows[source][role] = indices
            roles[role] = {
                "plates": role_plates,
                "profiles": len(indices),
                "membership_sha256": digest([[int(i), profiles[i]["source_row"],
                                               profiles[i]["profile_id"], profiles[i]["plate"]]
                                              for i in indices]),
            }
        partition = np.sort(np.concatenate(list(rows[source].values())))
        if not np.array_equal(partition, np.arange(len(profiles))):
            raise ValueError(f"{source}: split does not cover every profile exactly once")
        counts = tuple((len(roles[role]["plates"]), roles[role]["profiles"]) for role in ROLES)
        if not corpus.fixture and counts != EXPECTED_COUNTS[source]:
            raise ValueError(f"{source}: unexpected split counts {counts}")
        sources[source] = {"plate_order_sha256": digest(plates), "roles": roles}
    return {
        "schema_version": 1,
        "seed": SEED,
        "corpus_manifest_sha256": corpus.manifest_sha256,
        "gene_order_sha256": corpus.gene_order_sha256,
        "sources": sources,
    }, rows


def freeze_split(corpus, path):
    record, _ = split_memberships(corpus)
    publish_json(path, record)
    return record


class Split:
    def __init__(self, corpus, path):
        record_bytes = Path(path).read_bytes()
        record, self.rows = split_memberships(corpus)
        if json.loads(record_bytes) != record:
            raise ValueError("split does not match the frozen plate order, memberships or corpus")
        self.provenance = {
            "corpus_manifest_sha256": corpus.manifest_sha256,
            "gene_order_sha256": corpus.gene_order_sha256,
            "split_sha256": hashlib.sha256(record_bytes).hexdigest(),
            "normalization": corpus.normalization,
        }


def validate_batch_size(batch_size):
    if (not isinstance(batch_size, Integral) or isinstance(batch_size, bool)
            or not 1 <= batch_size <= MAX_BATCH_SIZE):
        raise ValueError(f"batch size must be between 1 and {MAX_BATCH_SIZE}")


def batches(corpus, split, source, role, batch_size):
    validate_batch_size(batch_size)
    rows = split.rows[source][role]
    for start in range(0, len(rows), batch_size):
        yield corpus.read(source, rows[start:start + batch_size])


def fit_mean(corpus, split, batch_size):
    total = np.zeros(corpus.gene_count, dtype=np.float64)
    count = 0
    for source in SOURCES:
        for values in batches(corpus, split, source, "train", batch_size):
            total += values.sum(axis=0, dtype=np.float64)
            count += len(values)
    return (total / count).astype(np.float32)


def artifact_payloads(directory):
    files = list(directory.iterdir())
    if any(not path.is_file() or path.is_symlink() for path in files):
        raise ValueError("artifact payloads must be regular files")
    return {path.name: file_sha256(path) for path in files if path.name != "manifest.json"}


@contextmanager
def publish_artifact(output, corpus, split, model_id, metadata):
    """Seal the selected inference files and publish them as one directory."""
    output = Path(output)
    if output.exists():
        raise ValueError(f"artifact already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=".artifact-") as pending:
        directory = Path(pending)
        yield directory
        payloads = artifact_payloads(directory)
        if not payloads or (directory / "manifest.json").exists():
            raise ValueError("artifact requires payload files and reserves manifest.json")
        record = {"schema_version": 1, "model_id": model_id, "gene_count": corpus.gene_count,
                  "provenance": split.provenance, "metadata": metadata, "files": payloads}
        (directory / "manifest.json").write_bytes(json_bytes(record))
        directory.rename(output)


def verify_artifact(path, corpus, split):
    directory = Path(path)
    manifest_path = directory / "manifest.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    if (manifest["schema_version"] != 1 or manifest["gene_count"] != corpus.gene_count
            or manifest["provenance"] != split.provenance):
        raise ValueError("artifact does not match the corpus, axis and frozen split")
    if not manifest["files"] or artifact_payloads(directory) != manifest["files"]:
        raise ValueError("artifact file digest mismatch")
    return {**manifest, "artifact_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest()}


def load_mean(path, corpus, split):
    record = verify_artifact(path, corpus, split)
    if record["model_id"] != "mean" or set(record["files"]) != {"mean.npy"}:
        raise ValueError("expected a training-gene mean artifact")
    mean = np.load(Path(path) / "mean.npy")
    if (mean.dtype != np.float32 or mean.shape != (corpus.gene_count,)
            or not np.isfinite(mean).all()):
        raise ValueError("invalid training-gene mean payload")
    return record, mean


def checked_test_role(corpus, role, allow_test):
    if role not in ROLES:
        raise ValueError(f"unknown split role: {role}")
    if role == "test" and not corpus.fixture and not allow_test:
        raise ValueError("real test scoring requires explicit authorization (--allow-test)")


def center(values):
    # Translating first gives exactly zero variance for a constant profile.
    centered = values - values[:, :1]
    return centered - centered.mean(axis=1, keepdims=True)


def score(corpus, split, role, predict, mean, batch_size, *, allow_test=False):
    checked_test_role(corpus, role, allow_test)
    mean = np.asarray(mean, dtype=np.float64)
    if mean.shape != (corpus.gene_count,) or not np.isfinite(mean).all():
        raise ValueError("invalid training-gene mean")
    results = {}
    for source in SOURCES:
        pcc_sum = rmse_sum = squared_error = reference_error = 0.0
        count = defined_count = 0
        for values in batches(corpus, split, source, role, batch_size):
            target = values.astype(np.float64)
            prediction = np.asarray(predict(values), dtype=np.float64)
            if prediction.shape != values.shape or not np.isfinite(prediction).all():
                raise ValueError(f"{source}: predictions must be finite and match the target axis and batch")
            if not np.isfinite(target).all():
                raise ValueError(f"{source}: nonfinite expression")
            sse = np.square(prediction - target).sum(axis=1)
            squared_error += float(sse.sum())
            reference_error += float(np.square(target - mean).sum(axis=1).sum())
            rmse_sum += float(np.sqrt(sse / corpus.gene_count).sum())
            observed, reconstructed = center(target), center(prediction)
            denominator = np.linalg.norm(observed, axis=1) * np.linalg.norm(reconstructed, axis=1)
            defined = denominator > 0
            dot = np.einsum("ij,ij->i", observed, reconstructed)
            pcc_sum += float(np.clip(dot[defined] / denominator[defined], -1, 1).sum())
            defined_count += int(defined.sum())
            count += len(values)
        results[source] = {
            "profiles": count,
            "defined_pcc_profiles": defined_count,
            "pcc": pcc_sum / defined_count if defined_count else None,
            "rmse": rmse_sum / count,
            "mse": squared_error / (count * corpus.gene_count),
            "relative_squared_error": squared_error / reference_error if reference_error else None,
        }
    return {"role": role, "sources": results}


def verify_replay(actual, expected):
    if not isinstance(expected, dict) or actual.keys() != expected.keys():
        raise ValueError("score replay fields differ")
    for key, observed in actual.items():
        previous = expected[key]
        if isinstance(observed, dict) and isinstance(previous, dict):
            verify_replay(observed, previous)
        elif isinstance(observed, float) and isinstance(previous, (int, float)):
            if not np.isclose(observed, previous, **REPLAY_TOLERANCE):
                raise ValueError(f"score replay differs for {key}: {observed} vs {previous}")
        elif observed != previous:
            raise ValueError(f"score replay differs for {key}: {observed} vs {previous}")
