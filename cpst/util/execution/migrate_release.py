"""Reissue accepted CPST metadata using the published original as the source."""

import argparse
from dataclasses import dataclass
from functools import lru_cache
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile
from tempfile import NamedTemporaryFile

STAGE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STAGE / "util/output/api"))

from cpst_util.release_names import current_path, migrate_names

ORIGINAL_REVISION = "9e7c85e38d2ab09c9b509ff8b2f30c67ea07ea85"
ORIGINAL_ROOT = "stage2_CarPS-VAE"
FIGURE_RECEIPTS = {"342-figure-render.json", "figure-render.json"}
CURRENT_STEPS = {"07_cpsv_alpha", "08_cpsv_chem-emb", "09_cpsv_axis", "10_cpsv_dose"}
LEDGERS = ("formal-receipts.sha256", "carps-vae-receipts.sha256")
ALPHA_MANIFEST = f"{ORIGINAL_ROOT}/07_cpsv_alpha/output/release.json"
FORMAL_CODE = {
    "08_cpsv_chem-emb": ("ac62a177a54aa99e6ac626728af12394f3e0d940", 54),
    "09_cpsv_axis": ("d3395ca0e0a055dc28e8e431740f5ec640ee4247", 63),
    "10_cpsv_dose": ("b8be1ad9244c9a3e03eba3e0c77ec64a59ed695d", 90),
}
METADATA_FIELDS = {".json", ".jsonl", ".csv", ".sha256"}
DIGEST = re.compile(rb"(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])")
RECORD = "cpst/util/execution/385-release-migration.json"
ARCHIVE = "cpst/util/output/migration-originals.tar.gz"


def original_blob(repo, name):
    return subprocess.check_output(
        ["git", "cat-file", "blob", f"{ORIGINAL_REVISION}:{name}"], cwd=repo,
    )


def digest(contents):
    return sha256(contents).hexdigest()


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def original_paths(repo):
    paths = subprocess.check_output(
        ["git", "ls-tree", "-rz", "--name-only", ORIGINAL_REVISION, "--", ORIGINAL_ROOT], cwd=repo,
    ).decode()
    return sorted(name for name in paths.split("\0") if name)


def migrate_figure(repo, name, original):
    receipt = json.loads(original)
    if receipt["schema"] != "stage2-figure-render-v1" or receipt["status"] != "PASS":
        raise ValueError(f"unaccepted figure receipt: {name}")
    original_root = str(Path(name).parent.parent) if name.endswith("/342-figure-render.json") else ORIGINAL_ROOT
    target_root = repo / current_path(original_root)
    renderers = {}
    for relative, expected in receipt["renderer_files_sha256"].items():
        if digest(original_blob(repo, f"{original_root}/{relative}")) != expected:
            raise ValueError(f"original figure source differs: {original_root}/{relative}")
        renamed = current_path(relative)
        renderers[renamed] = digest((target_root / renamed).read_bytes())
    receipt.update(
        schema="cpst-figure-render-v1",
        step=current_path(receipt["step"]),
        renderer_files_sha256=renderers,
    )
    return json_bytes(receipt)


@dataclass(frozen=True)
class Record:
    name: str
    accepted: bytes
    candidate: bytes
    tracked: bool
    dependencies: tuple[str, ...] = ()

    @property
    def original_hash(self):
        return digest(self.accepted)


def metadata_path(name):
    path = Path(name)
    if path.parts[:1] != (ORIGINAL_ROOT,) or ".." in path.parts or path.as_posix() != name:
        return False
    return (
        path.name in {*FIGURE_RECEIPTS, *LEDGERS}
        or (len(path.parts) > 3 and path.parts[1] in CURRENT_STEPS
            and ((path.parts[2] in {"log", "output"} and path.suffix in METADATA_FIELDS)
                 or (path.parts[2:4] == ("figure", "data") and path.suffix == ".csv")))
    )


def verify_original_ledgers(repo, read):
    """Verify the accepted pinned roots before computing any migrated digests."""
    for ledger in LEDGERS:
        contents = original_blob(repo, f"{ORIGINAL_ROOT}/util/training/{ledger}").decode()
        for row in contents.splitlines():
            expected, name = row.split()
            if digest(read(name)) != expected:
                raise ValueError(f"original acceptance ledger differs: {name}")


def verify_original_formal_code(repo, records):
    """Check each original formal registration's source dictionary at its pinned revision."""
    @lru_cache(maxsize=None)
    def blob_digest(revision, source):
        contents = subprocess.check_output(
            ["git", "cat-file", "blob", f"{revision}:{source}"], cwd=repo,
        )
        return digest(contents)

    for step, (revision, count) in FORMAL_CODE.items():
        prefix = f"{ORIGINAL_ROOT}/{step}/output/formal/matrix/runs/"
        registrations = [
            record for name, record in records.items()
            if name.startswith(prefix) and name.endswith("/registration.json")
        ]
        if len(registrations) != count:
            raise ValueError(f"expected {count} original formal registrations for {step}")
        for registration in registrations:
            original = json.loads(registration.accepted)["code_sha256"]
            if not original or any(blob_digest(revision, source) != expected for source, expected in original.items()):
                raise ValueError(f"historic source differs for {registration.name}")


def metadata_records(repo, tracked, names, read):
    records = {}
    for name in sorted(names):
        accepted = read(name)
        if Path(name).name in FIGURE_RECEIPTS:
            candidate = migrate_figure(repo, name, accepted)
        else:
            candidate = migrate_names(accepted)
            if name == ALPHA_MANIFEST:
                document = json.loads(candidate)
                # Step 06.4 retains its own original registration and exact Git paths.
                document["step_06_4"] = json.loads(accepted)["step_06_4"]
                candidate = json_bytes(document)
        records[name] = Record(name, accepted, candidate, name in tracked)
    return records


def dose_preflight(records):
    """Reissue the aggregate of the registered dose coordinates, in roster order."""
    matrix = f"{ORIGINAL_ROOT}/10_cpsv_dose/output/formal/matrix"
    name = f"{matrix}/preflight.json"
    record = records.get(name)
    if record is None:
        return
    reading = json.loads(record.accepted)
    names = tuple(f"{matrix}/runs/{arm}/fold-{fold}/seed-{seed}/registration.json"
                  for arm, fold, seed in reading["registered_roster"])

    def identity(registrations):
        configs = []
        for data in registrations:
            config = json.loads(data)
            del config["nonce"]
            configs.append(json.dumps(config, sort_keys=True))
        return digest(json.dumps(configs, separators=(",", ":")).encode())

    old_hash = identity(records[config].accepted for config in names)
    if reading["registered_configs_sha256"] != old_hash or reading["formal_configs_sha256"] != old_hash:
        raise ValueError("original Step 10 configuration aggregate differs")
    records[name] = Record(record.name, record.accepted, record.candidate, record.tracked, names)
    return identity


def build_plan(records):
    """Resolve metadata references in dependency order, with exact digest matches."""
    identity = dose_preflight(records)
    by_digest = {}
    for record in records.values():
        key = record.original_hash
        if key in by_digest:
            counterpart = records[by_digest[key]]
            if (record.candidate, record.dependencies) != (counterpart.candidate, counterpart.dependencies):
                raise ValueError(f"conflicting migration of identical original bytes: {record.name}")
        else:
            by_digest[key] = record.name
    resolved = {}
    visiting = set()

    def migrate(name):
        if name in resolved:
            return resolved[name]
        if name in visiting:
            raise ValueError(f"cycle in accepted metadata: {name}")
        visiting.add(name)
        record = records[name]

        def dependency(match):
            key = match[0].decode()
            return digest(migrate(by_digest[key])).encode() if key in by_digest else match[0]

        output = DIGEST.sub(dependency, record.candidate)
        if record.dependencies:
            current = json.loads(output)
            updated = identity(migrate(key) for key in record.dependencies)
            gate_hash = json.loads(record.accepted)["formal_configs_sha256"]
            current.update(
                registered_configs_sha256=updated, formal_configs_sha256=updated,
                gate_configs_sha256=gate_hash,
            )
            output = json_bytes(current)
        visiting.remove(name)
        resolved[name] = output
        return output

    for name in records:
        migrate(name)
    return resolved


def write_atomic(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(contents)
            handle.close()
            mode = path.stat().st_mode if path.exists() else 0o644
            temporary.chmod(mode)
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)


def migration_record(records, migrated):
    items = [
        {"original_path": name, "original_sha256": record.original_hash,
         "path": current_path(name), "sha256": digest(migrated[name]),
         "original_source": "git" if record.tracked else "sealed_archive"}
        for name, record in records.items()
        if migrated[name] != record.accepted or current_path(name) != name
    ]
    return {"schema": "cpst-release-migration-v1", "original_revision": ORIGINAL_REVISION, "metadata": items}


def write_archive(path, records):
    data = io.BytesIO()
    with gzip.GzipFile(fileobj=data, mode="wb", mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode="w") as archive:
            for name, record in records.items():
                if record.tracked:
                    continue
                entry = tarfile.TarInfo(name)
                entry.size = len(record.accepted)
                entry.mode = 0o644
                archive.addfile(entry, io.BytesIO(record.accepted))
    write_atomic(path, data.getvalue())


def verify(repo, tracked, archive_path, *, install=False):
    recorded = json.loads((repo / RECORD).read_bytes())
    if recorded["original_revision"] != ORIGINAL_REVISION:
        raise ValueError("unknown original CPST release revision")
    pending = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        def read(name):
            if name in tracked:
                return original_blob(repo, name)
            handle = archive.extractfile(name)
            if handle is None:
                raise ValueError(f"missing original metadata: {name}")
            with handle:
                return handle.read()

        verify_original_ledgers(repo, read)
        names = {entry["original_path"] for entry in recorded["metadata"]}
        if len(names) != len(recorded["metadata"]) or not all(map(metadata_path, names)):
            raise ValueError("invalid CPST metadata migration paths")
        records = metadata_records(repo, tracked, names, read)
        verify_original_formal_code(repo, records)
        migrated = build_plan(records)
        if recorded != migration_record(records, migrated):
            raise ValueError("CPST metadata digest correspondence differs from its sealed original")
        for name, expected in migrated.items():
            target = repo / current_path(name)
            observed = target.read_bytes()
            if observed == expected:
                continue
            if not install or name in tracked or observed != records[name].accepted:
                raise ValueError(f"migrated release metadata differs: {target}")
            if target.is_symlink() or not target.resolve().is_relative_to(repo):
                raise ValueError(f"migrated metadata cannot replace a linked payload: {target}")
            pending[target] = expected
    if install:
        if archive_path != repo / ARCHIVE:
            write_atomic(repo / ARCHIVE, archive_path.read_bytes())
        for target, contents in pending.items():
            write_atomic(target, contents)
    print(f"Verified {len(migrated)} CPST metadata files against their original release")


def migrate(repo, tracked, original):
    ignored = subprocess.check_output(
        ["git", "-C", str(original), "ls-files", "-zoi", "--exclude-standard", ORIGINAL_ROOT],
    ).decode().split("\0")

    def read(name):
        return original_blob(repo, name) if name in tracked else (original / name).read_bytes()

    verify_original_ledgers(repo, read)
    names = {name for name in tracked.union(ignored) if metadata_path(name)}
    records = metadata_records(repo, tracked, names, read)
    verify_original_formal_code(repo, records)
    migrated = build_plan(records)
    write_archive(repo / ARCHIVE, records)
    for name, expected in migrated.items():
        target = repo / current_path(name)
        if target.read_bytes() != expected:
            write_atomic(target, expected)
    record = migration_record(records, migrated)
    write_atomic(repo / RECORD, json_bytes(record))
    print(f"Reissued {len(record['metadata'])} CPST metadata files")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("migrate", "verify", "install"))
    parser.add_argument("--original-checkout", type=Path)
    parser.add_argument("--original-archive", type=Path)
    args = parser.parse_args()
    repo = STAGE.parent
    tracked = set(original_paths(repo))
    if args.action == "migrate":
        if args.original_checkout is None:
            parser.error("migrate requires --original-checkout")
        migrate(repo, tracked, args.original_checkout.resolve())
    else:
        archive = args.original_archive.resolve() if args.original_archive else repo / ARCHIVE
        verify(repo, tracked, archive, install=args.action == "install")


if __name__ == "__main__":
    main()
