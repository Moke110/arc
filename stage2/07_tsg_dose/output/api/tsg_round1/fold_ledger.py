"""Read the immutable Step 07 ADR 0078 member-level fold ledger."""
from dataclasses import dataclass
import gzip
from hashlib import file_digest, sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[3]
FOLD_ROOT = STEP / "output/frozen-folds"
FOLD_RECEIPT = STEP / "log/247-axis-local-smoke.json"
PAIR_SCHEMA = "tsg-dose-round1-axis-local-pair-v1"
SCENARIOS = ("cold-drug", "cold-drug-dose")
STRIP_PARTITIONS = {"excluded-held-parent", "excluded-held-dose"}
PARTITIONS = {"train", "validation", "test"} | STRIP_PARTITIONS


@dataclass(frozen=True)
class FoldLedger:
    partitions: pd.DataFrame
    identities: pd.DataFrame
    evidence: dict


def _pinned_json(path):
    payload = path.read_bytes()
    return json.loads(payload), sha256(payload).hexdigest()


def _check_doses(ledger):
    held_doses = []
    for fold in range(3):
        dose_fold = ledger.loc[ledger.scenario.eq("cold-drug-dose") & ledger.fold.eq(fold)]
        cold_fold = ledger.loc[ledger.scenario.eq("cold-drug") & ledger.fold.eq(fold)]
        cold_test = cold_fold.loc[cold_fold.partition.eq("test"), "drug_parent"]
        test_doses = dose_fold.loc[dose_fold.partition.eq("test"), "dose_um"].unique()
        if len(test_doses) != 1:
            raise ValueError(f"cold-drug-dose fold {fold} must contain one held dose")
        held_doses.append(float(test_doses[0]))
        held_parent = dose_fold.drug_parent.isin(cold_test)
        held_dose = dose_fold.dose_um.eq(test_doses[0])
        for label, expected in {
            "test": held_parent & held_dose,
            "excluded-held-parent": held_parent & ~held_dose,
            "excluded-held-dose": ~held_parent & held_dose,
        }.items():
            if not dose_fold.partition.eq(label).eq(expected).all():
                raise ValueError(f"incorrect {label} rows in cold-drug-dose fold {fold}")
    if len(set(held_doses)) != 3:
        raise ValueError("three distinct cold-drug-dose held doses are required")
    return held_doses


def load_folds(root, receipt_path, native_source, axis_sha256):
    """Verify the historical receipt chain and return ordered per-member labels."""
    root = Path(root)
    receipt, receipt_sha256 = _pinned_json(Path(receipt_path))
    pair, pair_sha256 = _pinned_json(root / "manifest.json")
    if (pair_sha256 != receipt["pair_sha256"] or pair["schema"] != PAIR_SCHEMA
            or pair["view_generation"]["axis"]["sha256"] != axis_sha256):
        raise ValueError("frozen Step 07 pair digest, schema or axis differs")
    expected_source = pair["sources"][native_source]
    source = root / native_source
    source_manifest, source_sha256 = _pinned_json(source / "manifest.json")
    if source_sha256 != expected_source:
        raise ValueError(f"frozen Step 07 source or fold digest differs: {native_source}")
    splits = source / "splits.csv"
    stream = splits.open("rb") if splits.is_file() else gzip.open(splits.with_suffix(".csv.gz"), "rb")
    with stream:
        splits_sha256 = file_digest(stream, "sha256").hexdigest()
        if splits_sha256 != source_manifest["outputs"]["splits.csv"]:
            raise ValueError(f"frozen Step 07 source or fold digest differs: {native_source}")
        stream.seek(0)
        ledger = pd.read_csv(stream, dtype={"profile_id": str, "drug_parent": str},
                             keep_default_na=False, float_precision="round_trip")
    expected_folds = [(scenario, fold) for scenario in SCENARIOS for fold in range(3)]
    if (list(ledger.columns) != ["row", "profile_id", "drug_parent", "dose_um",
                                "scenario", "fold", "partition"]
            or set(zip(ledger.scenario, ledger.fold)) != set(expected_folds)
            or not ledger.partition.isin(PARTITIONS).all()
            or ledger.loc[ledger.scenario.eq("cold-drug"), "partition"].isin(STRIP_PARTITIONS).any()
            or not pd.api.types.is_numeric_dtype(ledger.dose_um)
            or not np.isfinite(ledger.dose_um).all() or (ledger.dose_um <= 0).any()
            or ledger.duplicated(["profile_id", "scenario", "fold"]).any()
            or (ledger.profile_id.eq("") | ledger.drug_parent.eq("")).any()):
        raise ValueError(f"invalid frozen fold ledger: {native_source}")
    identities = ledger[["profile_id", "drug_parent", "dose_um"]].drop_duplicates().set_index("profile_id")
    if not identities.index.is_unique:
        raise ValueError("frozen member drug or dose differs across folds")
    partitions = ledger.pivot(index="profile_id", columns=["scenario", "fold"], values="partition")
    return FoldLedger(partitions.reindex(columns=expected_folds), identities, {
        "schema": PAIR_SCHEMA,
        "pair_sha256": pair_sha256,
        "source_manifest_sha256": source_sha256,
        "splits_sha256": splits_sha256,
        "receipt_sha256": receipt_sha256,
        "split_seed": receipt["input_source_evidence"][native_source]["split_seed"],
        "held_dose_um": _check_doses(ledger),
    })
