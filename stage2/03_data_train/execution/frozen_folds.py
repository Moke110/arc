"""Validate and carry forward the original Step 07 member-level split ledger."""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from stage2_util.derivation import release
from training_data import SCENARIOS, STRIPS, sha256

PAIR_SCHEMA = "tsg-dose-round1-axis-local-pair-v1"
PARTITIONS = {"train": "train", "validation": "valid", "test": "test",
              **dict.fromkeys(STRIPS, "excluded")}


@dataclass(frozen=True)
class FrozenFolds:
    partitions: pd.DataFrame
    identities: pd.DataFrame
    evidence: dict

    def group_labels(self, observations, groups):
        member_labels = self.partitions.reindex(observations.index)
        identities = self.identities.reindex(observations.index)
        if (member_labels.isna().any().any()
                or not identities.drug_parent.eq(observations.drug_parent_inchikey).all()
                or not identities.dose_um.eq(observations.dose_um).all()):
            raise ValueError("missing frozen fold assignment or mismatched member identity")
        per_group = member_labels.groupby(groups, sort=False)
        if not per_group.nunique().eq(1).all().all():
            raise ValueError("frozen fold assignments disagree between members of a pooled row")
        return [
            [PARTITIONS[label] for label in original] + [
                STRIPS.get(label, "") for label in original[3:]
            ]
            for original in per_group.first().itertuples(index=False, name=None)
        ]


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
    pair_path = root / "manifest.json"
    receipt = release.read_json(receipt_path)
    pair_sha256 = sha256(pair_path)
    pair = release.read_json(pair_path)
    if (pair_sha256 != receipt["pair_sha256"] or pair["schema"] != PAIR_SCHEMA
            or pair["view_generation"]["axis"]["sha256"] != axis_sha256):
        raise ValueError("frozen Step 07 pair digest, schema or axis differs")
    source = root / native_source
    source_path = source / "manifest.json"
    source_sha256 = sha256(source_path)
    source_manifest = release.read_json(source_path)
    splits = source / "splits.csv"
    splits_sha256 = sha256(splits)
    if (source_sha256 != pair["sources"][native_source]
            or splits_sha256 != source_manifest["outputs"]["splits.csv"]):
        raise ValueError(f"frozen Step 07 source or fold digest differs: {native_source}")
    expected_folds = [(scenario, fold) for scenario in SCENARIOS for fold in range(3)]
    ledger = pd.read_csv(splits, dtype={"profile_id": str, "drug_parent": str},
                         keep_default_na=False, float_precision="round_trip")
    if (list(ledger.columns) != ["row", "profile_id", "drug_parent", "dose_um",
                                "scenario", "fold", "partition"]
            or set(zip(ledger.scenario, ledger.fold)) != set(expected_folds)
            or not ledger.partition.isin(PARTITIONS).all()
            or ledger.loc[ledger.scenario.eq("cold-drug"), "partition"].isin(STRIPS).any()
            or not pd.api.types.is_numeric_dtype(ledger.dose_um)
            or not np.isfinite(ledger.dose_um).all() or (ledger.dose_um <= 0).any()
            or ledger.duplicated(["profile_id", "scenario", "fold"]).any()
            or (ledger.profile_id.eq("") | ledger.drug_parent.eq("")).any()):
        raise ValueError(f"invalid frozen fold ledger: {native_source}")
    identities = ledger[["profile_id", "drug_parent", "dose_um"]].drop_duplicates().set_index("profile_id")
    if not identities.index.is_unique:
        raise ValueError("frozen member drug or dose differs across folds")
    partitions = ledger.pivot(index="profile_id", columns=["scenario", "fold"], values="partition")
    return FrozenFolds(partitions.reindex(columns=expected_folds), identities, {
        "schema": PAIR_SCHEMA,
        "pair_sha256": pair_sha256,
        "source_manifest_sha256": source_sha256,
        "splits_sha256": splits_sha256,
        "receipt_sha256": sha256(receipt_path),
        "split_seed": receipt["input_source_evidence"][native_source]["split_seed"],
        "held_dose_um": _check_doses(ledger),
    })
