"""Read the versioned, shared tenv evaluation partitions."""

from collections import defaultdict
import csv
from dataclasses import dataclass
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path

from .seeds import read_seeds

COLD_DRUG = "tenv-cold-drug-v1"
COLD_DOSE_DRUG = "tenv-cold-dose-drug-v1"
SCHEMA = "cpst-tenv-evaluation-v1"
DEFAULT = Path(__file__).resolve().parents[2] / "tenv-folds" / "v1"
MANIFEST_SHA256 = "f6a6a9c43ec5e146fd3ee8a06c3aec62b98e686f18a0f0fd3920f6b4469b084c"
CONDITION_FILE = "conditions.csv"
CONDITION_COLUMNS = ("plate", "drug_parent", "dose_um")
HISTORICAL_COLUMNS = tuple(f"cold_drug_{index}" for index in range(3))
LEDGER_COLUMNS = (*CONDITION_COLUMNS, "parent_group", *HISTORICAL_COLUMNS)
PARTITIONS = ("train", "valid", "test")
HELD_DOSES = (0.2, 1.0, 10.0)
SPLIT_SEED = 20261128
COLD_SEEDS = (17, 42, 43, 2026, 2027)


@dataclass(frozen=True)
class Condition:
    plate: str
    drug_parent: str
    dose_um: float

    @classmethod
    def from_row(cls, row):
        return cls(row["plate"], row["drug_parent"], float(row["dose_um"]))


@dataclass(frozen=True)
class Fold:
    number: int
    parent_group: int
    held_parents: frozenset
    held_dose_um: float | None
    partitions: tuple

    def positions(self, conditions):
        """Return all split labels as positions in the complete caller row order."""
        labels = dict(self.partitions)
        conditions = tuple(conditions)
        if len(conditions) != len(labels) or set(conditions) != labels.keys():
            raise ValueError("tenv condition population differs from the frozen ledger")
        positions = defaultdict(list)
        for index, condition in enumerate(conditions):
            positions[labels[condition]].append(index)
        return {phase: tuple(rows) for phase, rows in positions.items()}

    def positions_for_treatments(self, treatments):
        """Align this fold to the ordered Step 03 treatment table."""
        return self.positions(Condition(*row) for row in
                              treatments[list(CONDITION_COLUMNS)].itertuples(index=False, name=None))


@dataclass(frozen=True)
class Protocol:
    name: str
    seeds: tuple
    folds: tuple
    release_sha256: str

    @property
    def coordinates(self):
        return tuple((fold.number, seed) for fold in self.folds for seed in self.seeds)


@dataclass(frozen=True)
class TenvFolds:
    release_sha256: str
    conditions: tuple
    axis_inputs: tuple
    protocols: tuple

    def protocol(self, name):
        """Select an immutable protocol by versioned ID."""
        for protocol in self.protocols:
            if protocol.name == name:
                return protocol
        raise ValueError(f"unknown tenv evaluation protocol: {name}")


def _parse_ledger(payload):
    reader = csv.DictReader(StringIO(payload.decode()))
    if tuple(reader.fieldnames or ()) != LEDGER_COLUMNS:
        raise ValueError("invalid tenv fold ledger columns")
    conditions = []
    parent_roles = {}
    for row in reader:
        condition = Condition.from_row(row)
        roles = (int(row["parent_group"]), *(row[name] for name in HISTORICAL_COLUMNS))
        if (not condition.plate or not condition.drug_parent or condition.dose_um not in HELD_DOSES
                or roles[0] not in range(5) or set(roles[1:]) - set(PARTITIONS)
                or roles[1:].count("test") != 1 or parent_roles.setdefault(condition.drug_parent, roles) != roles):
            raise ValueError("invalid tenv parent assignments")
        conditions.append(condition)
    if len(conditions) != 1383 or len(set(conditions)) != 1383 or len(parent_roles) != 459:
        raise ValueError("invalid tenv condition population")
    if sorted(sum(roles[0] == group for roles in parent_roles.values()) for group in range(5)) != [91, 92, 92, 92, 92]:
        raise ValueError("invalid five-group tenv parent balance")
    for group in range(3):
        counts = tuple(sum(roles[group + 1] == phase for roles in parent_roles.values()) for phase in PARTITIONS)
        if counts != (275, 31, 153):
            raise ValueError("invalid historical tenv parent groups")
    if len({(row.drug_parent, row.dose_um) for row in conditions}) != 459 * len(HELD_DOSES):
        raise ValueError("incomplete all-dose tenv condition population")
    return tuple(conditions), parent_roles


def _cold_drug_folds(conditions, roles):
    for group in range(5):
        validation = (group + 1) % 5
        labels = {parent: "test" if values[0] == group else "valid" if values[0] == validation else "train"
                  for parent, values in roles.items()}
        yield Fold(group, group, frozenset(parent for parent, role in labels.items() if role == "test"),
                   None, tuple((row, labels[row.drug_parent]) for row in conditions))


def _dose_label(role, row_dose, held_dose):
    if role == "test":
        return "test" if row_dose == held_dose else "excluded-held-parent"
    return "excluded-held-dose" if row_dose == held_dose else role


def _cold_dose_folds(conditions, roles):
    for group in range(3):
        held_parents = frozenset(parent for parent, values in roles.items() if values[group + 1] == "test")
        for dose_index, dose in enumerate(HELD_DOSES):
            yield Fold(group * len(HELD_DOSES) + dose_index, group, held_parents, dose, tuple(
                (row, _dose_label(roles[row.drug_parent][group + 1], row.dose_um, dose)) for row in conditions
            ))


def _build_protocols(conditions, roles, digest):
    return (
        Protocol(COLD_DRUG, COLD_SEEDS, tuple(_cold_drug_folds(conditions, roles)), digest),
        Protocol(COLD_DOSE_DRUG, read_seeds(), tuple(_cold_dose_folds(conditions, roles)), digest),
    )


def _partition_digests(protocols):
    digests = {}
    for protocol in protocols:
        payload = [
            [fold.number, fold.parent_group, fold.held_dose_um, sorted(fold.held_parents),
             [[row.plate, row.drug_parent, row.dose_um, label] for row, label in fold.partitions]]
            for fold in protocol.folds
        ]
        digests[protocol.name] = sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()
    return digests


def compute_partition_digests(ledger):
    """Materialize the canonical protocol memberships for a candidate ledger."""
    return _partition_digests(_build_protocols(*_parse_ledger(ledger), ""))


def open_tenv_folds(directory=None, *, manifest_sha256=None):
    """Verify the fold release and expose both protocols to numbered studies."""
    root = DEFAULT if directory is None else Path(directory)
    manifest_bytes = (root / "manifest.json").read_bytes()
    digest = sha256(manifest_bytes).hexdigest()
    expected = MANIFEST_SHA256 if manifest_sha256 is None else manifest_sha256
    if digest != expected:
        raise ValueError("tenv fold manifest digest differs")
    manifest = json.loads(manifest_bytes)
    protocols = {COLD_DRUG: {"folds": 5, "seeds": list(COLD_SEEDS)},
                 COLD_DOSE_DRUG: {"folds": 9, "seeds": list(read_seeds())}}
    if (manifest["schema"] != SCHEMA or manifest["protocols"] != protocols
            or manifest["split_seed"] != SPLIT_SEED or len(manifest["axis_inputs"]) != 8):
        raise ValueError("invalid tenv fold manifest")
    ledger = (root / CONDITION_FILE).read_bytes()
    if sha256(ledger).hexdigest() != manifest["outputs"][CONDITION_FILE]:
        raise ValueError("tenv condition ledger digest differs")
    conditions, roles = _parse_ledger(ledger)
    releases = _build_protocols(conditions, roles, digest)
    if _partition_digests(releases) != manifest["partition_sha256"]:
        raise ValueError("tenv protocol partition digests differ")
    return TenvFolds(digest, conditions, tuple(sorted(manifest["axis_inputs"].items())), releases)
