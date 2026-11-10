"""Build the Step 09 landmark predictor and evaluator releases from native counts."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
STEP = Path(__file__).resolve().parents[2]
STEP02 = REPO / "stage2/02_dataset_corpus"
STEP03 = REPO / "stage2/03_modes_tenv"
STEP04 = REPO / "stage2/04_pdiff_hyperparameter"
sys.path[:0] = [str(path / "execution/script") for path in (STEP03, STEP04)]
sys.path.insert(0, str(STEP02 / "execution"))

from cmdata.shared_cm import fingerprint, frozen_digests, gene_mapping, input_files, read_csv
from inputs import SEEDS, sha256, split_drugs, write_json
from tahoe_inputs import DRUGS as TAHOE_DRUGS, source_metadata

ROLES = ("train", "validation", "test")
SPLIT_SEED = 20261110
IDENTITY_SEED = 20261111
IDENTITY_WIDTH = 256
ROW_BATCH = 64
DEFAULT_NATIVE = STEP02 / "output/corpus/v2"
DEFAULT_TENV = STEP03 / "output/input"
DEFAULT_LANDMARKS = REPO / "pert_model/XPert/local/metadata/l1000_gene_info_978.csv"
DEFAULT_OUTPUT = STEP / "output/corpus/tahoe_tenv_bulk_l1000"


def file_hashes(root):
    return {str(path.relative_to(root)): sha256(path)
            for path in sorted(root.rglob("*")) if path.is_file()}


def landmark_axis(cm, author, crosswalk, fixture):
    author = author.assign(gene_idx=author.gene_idx.astype(int)).sort_values("gene_idx")
    if (not np.array_equal(author.gene_idx, np.arange(len(author)))
            or author.gene_id.duplicated().any()
            or (not fixture and (len(author) != 978 or len(cm) != 18164))):
        raise ValueError("invalid author landmark or CM axis")
    if any(cm[name].eq("").any() or cm[name].duplicated().any()
           for name in ("entrez_id", "hgnc_id", "symbol", "ensembl_id")):
        raise ValueError("ambiguous CM identity")
    extra = author.loc[~author.gene_id.isin(cm.entrez_id)]
    if extra.gene_id.tolist() != ["50810"]:
        raise ValueError("expected HDGFL3 as the only supplemental landmark")
    ensembl = extra.ensembl_id.item()
    hits = crosswalk.loc[crosswalk.source_gene_id.eq(ensembl), ["symbol", "hgnc_id"]]
    if len(hits.drop_duplicates()) != 1 or hits.symbol.iloc[0] != "HDGFL3":
        raise ValueError("ambiguous native HDGFL3 identity")
    supplement = pd.DataFrame([{"symbol": "HDGFL3", "entrez_id": "50810",
                               "hgnc_id": hits.hgnc_id.iloc[0], "ensembl_id": ensembl}])
    family = pd.concat([cm, supplement], ignore_index=True)
    family["identity_index"] = np.arange(len(family))
    return author[["gene_idx", "gene_id"]].merge(
        family, left_on="gene_id", right_on="entrez_id", validate="one_to_one"
    )


def tenv_cohort(obs, reference):
    drugs = pd.read_csv(reference / "drugs.csv", keep_default_na=False)
    folds = pd.read_csv(reference / "folds.csv", keep_default_na=False)
    pd.testing.assert_frame_equal(folds, split_drugs(drugs))
    if folds.groupby(["fold", "parent_connectivity"]).role.nunique().max() != 1:
        raise ValueError("tenv parent connectivity crosses roles")

    controls = pd.read_csv(reference / "controls.csv", keep_default_na=False)
    dmso = obs.is_control & ~obs.is_blank & obs.compound_name.str.upper().eq("DMSO") & obs.vehicle.str.upper().eq("DMSO")
    native_rows = pd.Series(np.arange(len(obs)), index=obs.index)
    if not np.array_equal(controls.control_index, np.arange(len(controls))):
        raise ValueError("tenv control order is not canonical")
    dmso_by_plate = {str(plate): native_rows.loc[part.index].to_numpy()
                     for plate, part in obs[dmso].groupby("plate", observed=True)}
    if set(dmso_by_plate) != set(controls.plate):
        raise ValueError("tenv DMSO plates differ from frozen control membership")
    control_rows = [dmso_by_plate[plate] for plate in controls.plate]
    if not np.array_equal(controls.n_wells, [len(rows) for rows in control_rows]):
        raise ValueError("tenv DMSO counts differ from frozen control membership")
    wells = pd.read_csv(reference / "wells.csv", keep_default_na=False)
    measured = obs[~obs.is_control & ~obs.is_blank]
    if wells.case_id.tolist() != measured.index.tolist():
        raise ValueError("tenv eligible well order differs from frozen membership")
    control_plate = controls.plate.to_numpy()[wells.control_index]
    if (not np.array_equal(wells.drug_key, measured.drug_key)
            or not np.array_equal(wells.dose_um, measured.dose_um)
            or not np.array_equal(control_plate, measured.plate)
            or not measured.time_h.eq(24.0).all()):
        raise ValueError("tenv well/dose/plate mismatch with frozen membership")
    wells["source_row"] = native_rows.loc[wells.case_id].to_numpy()
    wells["plate"] = control_plate
    wells = wells.merge(drugs[["drug_key", "parent_connectivity", "morgan_index"]],
                        on="drug_key", validate="many_to_one")
    cohorts = {}
    for fold in range(3):
        membership = folds.loc[folds.fold.eq(fold)].set_index("drug_key").role
        role_of_well = wells.drug_key.map(membership)
        if role_of_well.isna().any():
            raise ValueError(f"tenv fold {fold} omits a well")
        for role in ROLES:
            cases = wells.loc[role_of_well.eq(role)].reset_index(drop=True)
            reference_name = "test_cases" if role == "test" else role
            frozen = pd.read_csv(reference / f"fold_{fold}/{reference_name}.csv", keep_default_na=False)
            pd.testing.assert_frame_equal(cases[frozen.columns], frozen)
            if cases.empty:
                raise ValueError(f"empty tenv fold {fold} role {role}")
            cohorts[f"fold_{fold}/{role}"] = cases
    return cohorts, controls, control_rows, drugs, np.load(reference / "morgan.npy"), folds


def tahoe_cohort(obs, drug_metadata):
    cases, controls, control_rows, drugs, morgan, exclusions, counts = source_metadata(obs, drug_metadata)
    cases = cases.merge(drugs[["drug_key", "parent_connectivity", "morgan_index"]],
                        on="drug_key", validate="many_to_one")
    plates = sorted(cases.plate.unique(),
                    key=lambda plate: hashlib.sha256(f"{SPLIT_SEED}|{plate}".encode()).digest())
    if len(plates) < 2:
        raise ValueError("Tahoe requires at least two eligible plates")
    validation = set(plates[:max(1, len(plates) // 5)])
    mask = cases.plate.isin(validation)
    cohorts = {"train": cases[~mask], "validation": cases[mask]}
    return cohorts, controls, control_rows, drugs, morgan, exclusions, counts


def normalize_counts(counts):
    """Return post-projection log1p CP10K and the number of zero libraries."""
    values = counts.astype(np.float32)
    sums = counts.sum(axis=1, dtype=np.float64)
    scale = np.divide(10_000, sums, out=np.zeros_like(sums), where=sums != 0)
    values *= scale.astype(np.float32)[:, None]
    np.log1p(values, out=values)
    return values, int(np.count_nonzero(sums == 0))


def raw_batches(matrix, rows, columns):
    """Yield projected native integers in ascending, bounded native row windows."""
    if not len(rows) or np.any(rows[1:] <= rows[:-1]):
        raise ValueError("requested native row positions must be unique and increasing")
    sparse = isinstance(matrix, h5py.Group)
    if sparse:
        if matrix.attrs["encoding-type"] != "csr_matrix":
            raise ValueError("expected native CSR")
        ptr = matrix["indptr"][:]
        shape = tuple(matrix.attrs["shape"])
        data, indices = matrix["data"], matrix["indices"]
        if (ptr.shape != (shape[0] + 1,) or ptr[0] != 0 or ptr[-1] != len(data)
                or len(data) != len(indices) or np.any(np.diff(ptr) < 0)):
            raise ValueError("invalid CSR pointer table")
        dtype = data.dtype
    else:
        shape = matrix.shape
        dtype = matrix.dtype
    if dtype.kind not in "iu" or rows[0] < 0 or rows[-1] >= shape[0]:
        raise ValueError("invalid raw count matrix or requested rows")

    cuts = np.flatnonzero(np.diff(rows // ROW_BATCH)) + 1
    starts, stops = np.r_[0, cuts], np.r_[cuts, len(rows)]
    for low, high in zip(starts, stops, strict=True):
        requested = rows[low:high]
        start, stop = requested[0], requested[-1] + 1
        if sparse:
            lo, hi = ptr[start], ptr[stop]
            counts, col = data[lo:hi], indices[lo:hi]
            if np.any(col < 0) or np.any(col >= shape[1]):
                raise ValueError("invalid native CSR columns")
            slab = csr_matrix((counts, col, ptr[start:stop+1] - lo), shape=(stop-start, shape[1]))
            if not slab.has_canonical_format:
                raise ValueError("invalid native CSR column order or duplicates")
            projected = slab[requested-start][:, columns].toarray()
        else:
            counts = matrix[start:stop]
            projected = counts[requested-start][:, columns]
        if np.any(counts < 0):
            raise ValueError("negative native raw count")
        yield slice(low, high), projected


def write_expression(path, matrix, cases, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = cases.source_row.to_numpy()
    order = np.argsort(rows)
    mapped = np.lib.format.open_memmap(path, mode="w+", dtype=np.float32,
                                      shape=(len(rows), len(columns)))
    zero = 0
    for part, raw in raw_batches(matrix, rows[order], columns):
        values, zeros = normalize_counts(raw)
        mapped[order[part]] = values
        zero += zeros
    mapped.flush()
    return zero


def write_controls(path, matrix, control_rows, columns):
    rows = np.concatenate(control_rows)
    order = np.argsort(rows)
    indices = np.repeat(np.arange(len(control_rows)), [len(group) for group in control_rows])
    sums = np.zeros((len(control_rows), len(columns)), dtype=np.float64)
    zero = 0
    for part, raw in raw_batches(matrix, rows[order], columns):
        values, zeros = normalize_counts(raw)
        np.add.at(sums, indices[order[part]], values)
        zero += zeros
    np.save(path, (sums / np.bincount(indices)[:, None]).astype(np.float32))
    return zero


def write_axis_view(predictor, evaluator, source, matrix, mapping, cohorts, controls,
                    control_rows, drugs, morgan):
    destination = predictor / source
    destination.mkdir(parents=True)
    mapping.to_csv(destination / "mapping.csv", index=False)
    columns = mapping.native_column.to_numpy()
    if (morgan.shape != (len(drugs), 2048) or not np.isin(morgan, [0, 1]).all()
            or not np.array_equal(drugs.morgan_index, np.arange(len(drugs)))):
        raise ValueError(f"{source}: invalid Morgan fingerprint order or values")
    np.save(destination / "morgan.npy", morgan)
    controls.to_csv(destination / "controls.csv", index=False)
    # Keep exact matching well identities for evaluator-only raw-count contrasts.
    pd.DataFrame({"control_index": np.repeat(controls.control_index, [len(r) for r in control_rows]),
                  "source_row": np.concatenate(control_rows)}).to_csv(
                      evaluator / f"{source}_control_rows.csv", index=False)
    zeros = {"controls": write_controls(destination / "controls.npy", matrix, control_rows, columns)}
    test_zeros = {}
    for name, cases in cohorts.items():
        path = destination / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith("/test"):
            cases.drop(columns="source_row").to_csv(path.with_suffix(".csv"), index=False)
            path = evaluator / source / name
            path.parent.mkdir(parents=True, exist_ok=True)
            cases[["case_id", "source_row"]].to_csv(path.with_name("test_rows.csv"), index=False)
            test_zeros[name] = write_expression(path.with_name("test_expression.npy"), matrix, cases, columns)
        else:
            cases.to_csv(path.with_suffix(".csv"), index=False)
            zeros[name] = write_expression(path.with_name(f"{path.name}_expression.npy"),
                                           matrix, cases, columns)
    return zeros, test_zeros


def build(native, tenv_release, landmarks, tahoe_drugs, output, fixture=False):
    if output.exists():
        raise FileExistsError(output)
    files = {**input_files(native), "author_landmarks.csv": landmarks, "tahoe_drugs.parquet": tahoe_drugs}
    original = fingerprint(files)
    input_hashes = {name: sha256(path) for name, path in files.items()}
    frozen_paths = ["controls.csv", "drugs.csv", "folds.csv", "morgan.npy", "wells.csv",
                    *(f"fold_{fold}/{name}.csv" for fold in range(3)
                      for name in ("train", "validation", "test_cases"))]
    tenv_hashes = {name: sha256(tenv_release / name) for name in frozen_paths}
    if not fixture:
        for name, digest in frozen_digests().items():
            if input_hashes[name] != digest:
                raise ValueError(f"native input differs from pinned Step 02 digest: {name}")
        tenv_manifest = tenv_release / "manifest.json"
        pin = json.loads((STEP03 / "execution/input_record/verification.json").read_text())["manifest_sha256"]
        if sha256(tenv_manifest) != pin:
            raise ValueError("tenv release differs from the pinned Step 03 manifest")
        expected = json.loads(tenv_manifest.read_text())["files"]
        if any(expected[name] != digest for name, digest in tenv_hashes.items()):
            raise ValueError("tenv release payload differs from pinned Step 03 artifacts")
    axis, crosswalk = read_csv(files["genes.csv"]), read_csv(files["gene_crosswalk.csv"])
    ordered = landmark_axis(axis, read_csv(landmarks), crosswalk, fixture)
    rng = np.random.default_rng(IDENTITY_SEED)
    vectors = rng.standard_normal((len(axis) + 1, IDENTITY_WIDTH), dtype=np.float32)
    if len(np.unique(vectors, axis=0)) != len(vectors):
        raise ValueError("duplicate frozen gene identity vector")

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=f".{output.name}-") as temporary:
        root = Path(temporary) / "release"
        predictor, evaluator = root / "predictor", root / "evaluator"
        predictor.mkdir(parents=True)
        evaluator.mkdir()
        shutil.copyfile(files["genes.csv"], predictor / "cm_genes.csv")
        ordered.to_csv(predictor / "genes.csv", index=False)
        np.save(predictor / "identity.npy", vectors[:-1])
        np.save(predictor / "identity_hdgfl3.npy", vectors[-1])
        zeros, test_zeros, cohorts_record = {}, {}, {}
        for source in ("tahoe", "tenv"):
            with h5py.File(files[f"{source}.h5ad"], "r") as handle:
                obs, var = read_elem(handle["obs"]), read_elem(handle["var"])
                if obs.index.has_duplicates:
                    raise ValueError(f"duplicate {source} observation identities")
                mapping = gene_mapping(source, ordered, var, crosswalk, fixture=True)
                if source == "tenv":
                    cohorts, controls, rows, drugs, morgan, folds = tenv_cohort(obs, tenv_release)
                    folds.to_csv(predictor / "tenv_folds.csv", index=False)
                else:
                    cohorts, controls, rows, drugs, morgan, exclusions, counts = tahoe_cohort(
                        obs, pd.read_parquet(tahoe_drugs)
                    )
                    exclusions.to_csv(predictor / "tahoe_exclusions.csv", index=False)
                    write_json(predictor / "tahoe_population.json", counts)
                zeros[source], test_zeros[source] = write_axis_view(
                    predictor, evaluator, source, handle["X"], mapping, cohorts, controls,
                    rows, drugs, morgan
                )
                cohorts_record[source] = {name: {"wells": len(cases), "groups": cases.group_id.nunique(),
                                                 "plates": sorted(cases.plate.unique())}
                                          for name, cases in cohorts.items()}
        if original != fingerprint(files):
            raise ValueError("native input changed during the corpus build")
        write_json(predictor / "manifest.json", {
            "schema": "step09-predictor-v1", "axis": "landmark", "gene_count": len(ordered),
            "fixture": fixture, "identity_seed": IDENTITY_SEED, "identity_generator": "NumPy PCG64",
            "identity_width": IDENTITY_WIDTH, "tahoe_split_seed": SPLIT_SEED, "training_seeds": SEEDS,
            "source_hashes": input_hashes, "step03_payload_sha256": tenv_hashes,
            "builder_sha256": sha256(Path(__file__)), "cohorts": cohorts_record,
            "zero_libraries": zeros, "files": file_hashes(predictor),
        })
        write_json(evaluator / "manifest.json", {
            "schema": "step09-evaluator-v1", "predictor_manifest_sha256": sha256(predictor / "manifest.json"),
            "test_zero_libraries": test_zeros, "files": file_hashes(evaluator),
        })
        root.rename(output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, default=DEFAULT_NATIVE)
    parser.add_argument("--tenv-release", type=Path, default=DEFAULT_TENV)
    parser.add_argument("--landmarks", type=Path, default=DEFAULT_LANDMARKS)
    parser.add_argument("--tahoe-drugs", type=Path, default=REPO / TAHOE_DRUGS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--fixture", action="store_true", help="allow compact, unpinned test inputs")
    args = parser.parse_args()
    print(build(args.native, args.tenv_release, args.landmarks, args.tahoe_drugs, args.output, args.fixture))


if __name__ == "__main__":
    main()
