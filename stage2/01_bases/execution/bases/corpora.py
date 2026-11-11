"""Bounded build and raw-count replay of the two indexed training corpora.

The training reader opens only one h5ad; it never discovers or opens truth labels.
Legacy shared expression is verified before reuse, never relabelled in place.
"""

from collections import defaultdict
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

BATCH = 64
NORMALIZATION = "project raw counts to frozen axis; CP10K; natural log1p; mean normalized matched controls"
SCHEMA = "stage2-indexed-training-corpora-v1"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def shared_library():
    root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    sys.path.insert(0, str(root / "stage2/02_dataset_corpus/execution"))
    from cmdata import shared_cm

    return shared_cm


@lru_cache(maxsize=4096)
def parent_structure(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError("invalid_structure")
    exact = Chem.MolToInchiKey(mol)
    parent = rdMolStandardize.Uncharger().uncharge(
        rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol))
    )
    Chem.AssignStereochemistry(parent, cleanIt=True, force=True)
    key = Chem.MolToInchiKey(parent)
    if len(key) != 27:
        raise ValueError("unresolved_parent_structure")
    return key, Chem.MolToSmiles(parent, isomericSmiles=True), exact


def dose_micromolar(value, unit="uM"):
    factors = {"uM": 1.0, "µM": 1.0, "μM": 1.0, "nM": 0.001, "mM": 1000.0, "M": 1e6}
    if unit not in factors:
        raise ValueError("unresolvable_dose_unit")
    dose = float(value) * factors[unit]
    if not np.isfinite(dose) or dose <= 0:
        raise ValueError("invalid_treated_dose")
    return dose


def metadata(source, obs, profiles, expression, drugs):
    """Resolve source rows, exact matched controls and administered parent evidence."""
    if obs.index.has_duplicates or profiles.source_row.duplicated().any():
        raise ValueError("duplicate source profiles")
    profile_position = dict(zip(profiles.source_row, range(len(profiles)), strict=True))
    structures = {}
    if source == "tahoe":
        for name, frame in drugs.dropna(subset=["drug"]).groupby("drug", sort=True):
            values = sorted(
                set(
                    str(x).strip()
                    for x in frame.canonical_smiles.dropna()
                    if str(x).strip()
                )
            )
            resolved = []
            for value in values:
                try:
                    resolved.append((value, *parent_structure(value)))
                except ValueError:
                    pass
            if len({x[1] for x in resolved}) > 1:
                raise ValueError(f"ambiguous Tahoe parent: {name}")
            structures[str(name).strip()] = resolved
    if source == "tenv":
        conditions = defaultdict(set)
        for row in obs.to_dict("records"):
            if row["is_control"] or row.get("is_blank", False):
                continue
            try:
                parent, _, _ = parent_structure(row["smiles"])
                dose = dose_micromolar(row["dose_um"])
            except (ValueError, TypeError):
                continue
            conditions[parent, dose].add(
                (str(row.get("time_h", "")), str(row.get("vehicle", "")))
            )
        conflicts = {
            f"{key[0]}|{key[1]:.12g}": sorted(values)
            for key, values in conditions.items()
            if len(values) != 1
            or any(x in ("", "nan", "None") for value in values for x in value)
        }
        if conflicts:
            raise ValueError(
                "incompatible or unresolved tenv parent-dose conditions: "
                + json.dumps(conflicts)
            )
    controls = defaultdict(list)
    for row in profiles.itertuples(index=False):
        native = obs.iloc[int(row.source_row)]
        control = bool(native.is_control) and str(
            native["drug" if source == "tahoe" else "compound_name"]
        ).strip() == ("DMSO_TF" if source == "tahoe" else "DMSO")
        if control:
            key = (
                (str(native.plate), str(native.cell_line_id))
                if source == "tahoe"
                else (str(native.plate),)
            )
            values = expression[profile_position[row.source_row]]
            if np.isfinite(values).all() and (values >= 0).all() and np.any(values > 0):
                controls[key].append(int(row.source_row))
    samples, exclusions = [], []
    for source_row, native in enumerate(obs.to_dict("records")):
        identifier = str(obs.index[source_row])
        if bool(native.get("is_blank", False)):
            exclusions.append(
                dict(
                    source_row=source_row,
                    source_profile_id=identifier,
                    reason="technical_blank",
                )
            )
            continue
        if bool(native["is_control"]):
            continue
        reason = None
        try:
            dose = dose_micromolar(native["dose_um"])
            if source == "tahoe":
                candidates = structures.get(str(native["drug"]).strip(), [])
                if not candidates:
                    raise ValueError("missing_or_invalid_structure")
                smiles, parent, parent_smiles, exact = candidates[0]
            else:
                smiles = native["smiles"]
                parent, parent_smiles, exact = parent_structure(smiles)
            key = (
                (str(native["plate"]), str(native["cell_line_id"]))
                if source == "tahoe"
                else (str(native["plate"]),)
            )
            matched = controls.get(key, [])
            if not matched:
                raise ValueError("missing_matched_control")
            if source == "tenv":
                for field in ("time_h", "vehicle"):
                    values = {str(obs.iloc[row][field]) for row in matched}
                    if values != {str(native[field])}:
                        raise ValueError("incompatible_matched_control_conditions")
            position = profile_position.get(source_row)
            if position is None:
                raise ValueError("missing_expression_profile")
            values = expression[position]
            if (
                not np.isfinite(values).all()
                or (values < 0).any()
                or not np.any(values > 0)
            ):
                raise ValueError("invalid_projected_profile")
        except (ValueError, TypeError, KeyError) as error:
            reason = str(error)
        if reason:
            exclusions.append(
                dict(source_row=source_row, source_profile_id=identifier, reason=reason)
            )
            continue
        record = {
            k: v for k, v in native.items() if isinstance(v, (str, int, float, bool))
        }
        record.update(
            sample_id=f"{source}:{identifier}",
            source=source,
            source_row=source_row,
            source_profile_id=identifier,
            expression_row=position,
            drug_parent=parent,
            parent_smiles=parent_smiles,
            source_smiles=str(smiles),
            source_inchikey=exact,
            source_parent_key=str(native.get("drug_key", "")),
            parent_identity_method="RDKit Cleanup; FragmentParent; Uncharger; stereochemistry retained",
            dose_um=dose,
            dose_original_value=float(native["dose_um"]),
            dose_original_unit="uM",
            dose_conversion="native dose_um * 1; upstream dose_source retained",
            cell_identity="hiPSC-CM"
            if source == "tenv"
            else str(native.get("cell_name", native["cell_line_id"])),
            matched_control_rows=json.dumps(matched, separators=(",", ":")),
            matched_control_ids=json.dumps(
                [str(obs.index[x]) for x in matched], separators=(",", ":")
            ),
            group_id=hashlib.sha256(f"{parent}|{dose:.12g}".encode()).hexdigest()[:24],
        )
        samples.append(record)
    if not samples:
        raise ValueError(f"{source}: no eligible samples")
    frame = pd.DataFrame(samples).set_index("sample_id", drop=False)
    frame.index.name = "sample_index"
    if frame.index.has_duplicates:
        raise ValueError("duplicate sample_id")
    if source == "tenv":
        conflicts = []
        for key, part in frame.groupby(["drug_parent", "dose_um"], sort=True):
            for field in ("time_h", "vehicle"):
                if (
                    field not in part
                    or part[field].isna().any()
                    or part[field].astype(str).eq("").any()
                    or part[field].nunique() != 1
                ):
                    conflicts.append(
                        {
                            "parent_dose": key,
                            "field": field,
                            "samples": part.sample_id.tolist(),
                        }
                    )
        if conflicts:
            raise ValueError(
                "incompatible or unresolved tenv parent-dose conditions: "
                + json.dumps(conflicts)
            )
    excluded = pd.DataFrame(
        exclusions, columns=["source_row", "source_profile_id", "reason"]
    )
    summary = {
        "samples": len(frame),
        "raw_profiles": len(obs),
        "declared_treated": int((~obs.is_control & ~obs.get("is_blank", False)).sum()),
        "controls": int(obs.is_control.sum()),
        "matched_control_profiles": len(set(x for v in controls.values() for x in v)),
        "excluded_by_reason": {
            str(k): int(v)
            for k, v in excluded.reason.value_counts().sort_index().items()
        },
        "parent_identity_remapped_samples": int(
            frame.source_parent_key.ne("")
            .mul(frame.source_parent_key.ne(frame.drug_parent))
            .sum()
        ),
        "parent_dose_groups": int(
            frame[["drug_parent", "dose_um"]].drop_duplicates().shape[0]
        ),
    }
    return frame, excluded, controls, summary


def means(expression, profiles, controls):
    position = dict(zip(profiles.source_row, range(len(profiles)), strict=True))
    return {
        key: np.mean(
            expression[[position[x] for x in rows]], axis=0, dtype=np.float64
        ).astype("float32")
        for key, rows in controls.items()
    }


def write_h5ad(path, frame, axis, expression, control_means, source):
    with h5py.File(path, "w") as output:
        output.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        write_elem(output, "obs", frame)
        var = axis.copy()
        var.insert(0, "gene_index", np.arange(len(axis), dtype="int32"))
        var.index = pd.Index(np.arange(len(axis)).astype(str), name="gene_index_label")
        write_elem(output, "var", var)
        write_elem(
            output,
            "uns",
            {"schema": SCHEMA, "normalization": NORMALIZATION, "source": source},
        )
        layers = output.create_group("layers")
        layers.attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
        arrays = {}
        for key, parent in [("X", output), ("control", layers), ("delta", layers)]:
            data = parent.create_dataset(
                key,
                shape=(len(frame), len(axis)),
                dtype="float32",
                chunks=(min(BATCH, len(frame)), len(axis)),
                compression="gzip",
                compression_opts=1,
                shuffle=True,
            )
            data.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})
            arrays[key] = data
        for start in range(0, len(frame), BATCH):
            part = frame.iloc[start : start + BATCH]
            treated = np.asarray(
                expression[part.expression_row.to_numpy()], dtype="float32"
            )
            control = np.stack(
                [
                    control_means[
                        (str(row.plate), str(row.cell_line_id))
                        if source == "tahoe"
                        else (str(row.plate),)
                    ]
                    for row in part.itertuples()
                ]
            )
            arrays["X"][start : start + len(part)] = treated
            arrays["control"][start : start + len(part)] = control
            arrays["delta"][start : start + len(part)] = treated - control


def training_header(path):
    """Read small metadata only; this API cannot access a truth-label release."""
    with h5py.File(path, "r") as file:
        if read_elem(file["uns"])["schema"] != SCHEMA:
            raise ValueError("unsupported training corpus schema")
        if set(file["layers"]) != {"control", "delta"}:
            raise ValueError("unexpected training layers")
        return read_elem(file["obs"]), read_elem(file["var"])


def replay_source(native, release, source, frame, controls, axis, crosswalk, fixture):
    """Recompute every treated/control value from native raw counts, not shared expression."""
    library = shared_library()
    # Only matched-control means are retained (<=700 strata); raw rows stream once.
    sums = {key: np.zeros(len(axis), dtype="float64") for key in controls}
    control_key = {row: key for key, rows in controls.items() for row in rows}
    counts_per = defaultdict(int)
    with h5py.File(native / f"{source}.h5ad") as raw:
        profiles, mapping, shape, _ = library.source_contract(
            raw, source, axis, crosswalk, fixture
        )
        rows = profiles.source_row.to_numpy()
        columns = mapping.native_column.to_numpy()
        for offset, counts, totals in library.projected_batches(
            raw["X"], shape, rows, columns, source
        ):
            expected = np.log1p(counts * (10000 / totals[:, None])).astype("float32")
            for index, row in enumerate(rows[offset : offset + len(counts)]):
                if row in control_key:
                    key = control_key[row]
                    sums[key] += expected[index]
                    counts_per[key] += 1
        control_means = {
            k: (v / counts_per[k]).astype("float32") for k, v in sums.items()
        }
        target_positions = {int(r): i for i, r in enumerate(frame.source_row)}
        max_error = 0.0
        with h5py.File(release, "r") as output:
            actual_obs, actual_var = training_header(release)
            pd.testing.assert_frame_equal(
                actual_obs, frame, check_dtype=False, check_categorical=False
            )
            if list(actual_var.gene_index) != list(range(len(axis))) or list(
                actual_var.hgnc_id
            ) != list(axis.hgnc_id):
                raise ValueError("training gene axis differs")
            if output["X"].shape != (len(frame), len(axis)):
                raise ValueError("training expression shape differs")
            replayed = 0
            for offset, counts, totals in library.projected_batches(
                raw["X"], shape, rows, columns, source
            ):
                expected = np.log1p(counts * (10000 / totals[:, None])).astype(
                    "float32"
                )
                selected = [
                    (local, target_positions[int(row)])
                    for local, row in enumerate(rows[offset : offset + len(counts)])
                    if int(row) in target_positions
                ]
                if not selected:
                    continue
                local, positions = map(np.array, zip(*selected, strict=True))
                part = frame.iloc[positions]
                control = np.stack(
                    [
                        control_means[
                            (str(row.plate), str(row.cell_line_id))
                            if source == "tahoe"
                            else (str(row.plate),)
                        ]
                        for row in part.itertuples()
                    ]
                )
                for name, reference in [
                    ("X", expected[local]),
                    ("layers/control", control),
                    ("layers/delta", expected[local] - control),
                ]:
                    actual = output[name][positions[0] : positions[-1] + 1]
                    if len(actual) != len(positions):
                        raise ValueError("source order is not contiguous")
                    max_error = max(
                        max_error, float(np.max(np.abs(actual - reference)))
                    )
                    if not np.allclose(actual, reference, atol=2e-6, rtol=1e-6):
                        raise ValueError(
                            f"{source} raw source replay differs: {name} row{positions[0]}"
                        )
                replayed += len(positions)
    return {
        "rows_replayed": replayed,
        "max_abs_error": max_error,
        "control_profiles_replayed": len(control_key),
    }


def publish(native, shared, drug_metadata, output, *, fixture=False, verify=False):
    library = shared_library()
    if fixture and not shared.exists():
        library.build(native, shared, True)
    # Real native bytes and every shared payload are pinned before reuse.
    shared_check = library.verify(native, shared, fixture)
    axis = library.read_csv(native / "tables/genes.csv")
    crosswalk = library.read_csv(native / "tables/gene_crosswalk.csv")
    drugs = pd.read_parquet(drug_metadata)
    prepared = {}
    summaries = {}
    for source in ("tahoe", "tenv"):
        with h5py.File(native / f"{source}.h5ad", "r") as file:
            obs = read_elem(file["obs"])
        profiles = library.read_csv(shared / f"{source}_profiles.csv")
        profiles["source_row"] = profiles.source_row.astype(int)
        expression = np.load(shared / f"{source}_expression.npy", mmap_mode="r")
        frame, excluded, controls, summary = metadata(
            source, obs, profiles, expression, drugs
        )
        prepared[source] = (frame, excluded, controls, profiles, expression)
        summaries[source] = summary
    inputs = {
        "shared_manifest": sha256(shared / "manifest.json"),
        "drug_metadata": sha256(drug_metadata),
    }
    manifest = {
        "schema": SCHEMA,
        "fixture": fixture,
        "axis_size": len(axis),
        "normalization": NORMALIZATION,
        "inputs": inputs,
        "rdkit_version": rdBase.rdkitVersion,
        "sources": summaries,
        "historical_inputs": "unchanged; native sources and prior experiments retain original semantics",
    }
    if verify:
        stored = json.loads((output / "manifest.json").read_text())
        for key, value in manifest.items():
            if stored[key] != value:
                raise ValueError(f"release metadata differs: {key}")
        for name, digest in stored["files"].items():
            if sha256(output / name) != digest:
                raise ValueError(f"release hash differs: {name}")
        checks = {}
        for source, (
            frame,
            excluded,
            controls,
            profiles,
            expression,
        ) in prepared.items():
            name = "tahoe_bulk" if source == "tahoe" else "tenv"
            if (output / f"{name}_exclusions.csv").read_text() != excluded.to_csv(
                index=False
            ):
                raise ValueError("exclusion accounting differs")
            checks[source] = replay_source(
                native,
                output / f"{name}.h5ad",
                source,
                frame,
                controls,
                axis,
                crosswalk,
                fixture,
            )
        return {
            "result": "PASS",
            "manifest_sha256": sha256(output / "manifest.json"),
            "sources": checks,
            "shared_source_replay": shared_check,
        }
    if output.exists():
        raise ValueError(f"release exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".training-", dir=output.parent) as tmp:
        stage = Path(tmp)
        names = []
        for source, (
            frame,
            excluded,
            controls,
            profiles,
            expression,
        ) in prepared.items():
            name = "tahoe_bulk" if source == "tahoe" else "tenv"
            write_h5ad(
                stage / f"{name}.h5ad",
                frame,
                axis,
                expression,
                means(expression, profiles, controls),
                source,
            )
            excluded.to_csv(stage / f"{name}_exclusions.csv", index=False)
            names.extend([f"{name}.h5ad", f"{name}_exclusions.csv"])
        manifest["files"] = {name: sha256(stage / name) for name in names}
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        os.rename(stage, output)
    return manifest
