"""Locked phase transitions and atomic evidence for registered study runs."""

from contextlib import contextmanager
from dataclasses import dataclass
import os
from pathlib import Path
import resource
from tempfile import TemporaryDirectory
import time

import torch

from stage2_util.derivation.release import read_json, sha256, write_json
from stage2_util.publication import publication_lock
from .transigen import PROCESS, PHASES, file_hashes, verified_receipt


def register(output, bind):
    with publication_lock(output) as root:
        if root.exists():
            raise ValueError("run registration requires an unused output directory")
        config = bind()
        root.mkdir()
        path = root / "registration.json"
        write_json(path, {**config, "nonce": PROCESS})
        return sha256(path)


@dataclass(frozen=True)
class Phase:
    name: str
    root: Path
    scratch: Path
    registration_sha: str
    config: dict
    training: dict | None
    sealed: dict | None
    files: tuple
    start: float

    def publish(self, details, identity, *, report_replay_resources=False):
        if self.sealed is not None and not report_replay_resources:
            return self.sealed
        if self.name == "replay":
            details = {**details, "test_replay": "exact"}
            hashes = file_hashes(self.root, self.files)
        else:
            hashes = file_hashes(self.scratch, self.files)
            for name in self.files:
                (self.scratch / name).replace(self.root / name)
        record = {
            **details, **identity, "outputs": hashes, "registration_sha256": self.registration_sha,
            "process_token": PROCESS, "pid": os.getpid(), "wall_seconds": round(time.monotonic() - self.start, 3),
            "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "cuda_peak_allocated_bytes": torch.cuda.max_memory_allocated(),
        }
        if self.sealed is None:
            name = PHASES[self.name]
            write_json(self.scratch / name, record)
            (self.scratch / name).replace(self.root / name)
        return record


@contextmanager
def phase(output, name, files, *, separate_registration=False):
    start = time.monotonic()
    with publication_lock(output) as root:
        registration = root / "registration.json"
        config = read_json(registration)
        registration_sha = sha256(registration)
        if separate_registration and config["nonce"] == PROCESS:
            raise ValueError("run phase must be separate from registration")
        receipt = PHASES[name]
        exists = (root / receipt).exists()
        if exists and name != "replay":
            raise ValueError(f"run phase already recorded: {name}")
        sealed = verified_receipt(root, receipt, registration_sha, files[name]) if exists else None
        training = verified_receipt(root, PHASES["train"], registration_sha, files["train"]) if name != "train" else None
        if name == "replay":
            verified_receipt(root, PHASES["predict"], registration_sha, files["predict"])
        with TemporaryDirectory(prefix=".phase-", dir=root) as temporary:
            yield Phase(name, root, Path(temporary), registration_sha, config, training, sealed, files[name], start)
