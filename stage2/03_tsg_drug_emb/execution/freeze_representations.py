"""Extract and verify the full seven-arm tenv parent representation freeze."""
import argparse
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / "output/api"))
from tsg_drug_emb.representations import (
    ARMS, PARENT_COLUMNS, SCHEMA, WIDTHS, RepresentationInputs, checked_vectors,
    fold_coverage, sha256,
)

KPGT_BASE_SHA256 = "8811c117cfedddff5ef2a512290d1a7a73d3d9560773860bc3b06942c597d15d"
CIRCULAR = {"ECFP4C": (2, False), "ECFP6C": (3, False),
            "FCFP4C": (2, True), "FCFP6C": (3, True)}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def circular(parents, directory):
    molecules = [Chem.MolFromSmiles(smiles) for smiles in parents.drug_parent_smiles]
    provenance = {}
    for arm, (radius, feature_atoms) in CIRCULAR.items():
        kwargs = {"atomInvariantsGenerator": rdFingerprintGenerator.GetMorganFeatureAtomInvGen()} if feature_atoms else {}
        generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, includeChirality=True, fpSize=WIDTHS[arm], **kwargs)
        matrix = np.asarray([generator.GetFingerprintAsNumPy(mol) for mol in molecules],
                            dtype=np.float32)
        np.save(directory / f"{arm}.npy", matrix)
        provenance[arm] = dict(radius=radius, feature_atoms=feature_atoms, include_chirality=True,
                               rdkit=rdBase.rdkitVersion)
    return provenance


def e3fp(parents, directory):
    from e3fp.pipeline import fprints_from_smiles

    matrix = np.full((len(parents), WIDTHS["E3FP"]), np.nan, dtype=np.float32)
    counts = np.zeros(len(parents), dtype=int)
    failures = {}
    conformer = dict(seed=42, first=3)
    fingerprint = dict(bits=WIDTHS["E3FP"], first=3, counts=False)
    for i, (key, smiles) in enumerate(parents.itertuples(index=False, name=None)):
        try:
            fps = fprints_from_smiles(smiles, key, confgen_params=conformer,
                                      fprint_params=fingerprint)
            if not 1 <= len(fps) <= 3:
                raise ValueError(f"unexpected E3FP conformer count {len(fps)}")
            bits = np.zeros(WIDTHS["E3FP"], dtype=bool)
            for fp in fps:
                bits[fp.indices] = True
            matrix[i] = bits
            counts[i] = len(fps)
        except (RuntimeError, ValueError, AttributeError) as error:
            failures[key] = str(error)
    np.save(directory / "E3FP.npy", matrix)
    parents[["drug_parent_inchikey"]].assign(conformers=counts).to_csv(
        directory / "e3fp_conformers.csv", index=False)
    return dict(e3fp=version("e3fp"), rdkit=rdBase.rdkitVersion, conformer=conformer,
                fingerprint=fingerprint, aggregation="bitwise OR", failures=failures)


def unimol(parents, directory, args):
    import torch
    from unimol_tools import UniMolRepr

    torch.set_num_threads(2)
    model = UniMolRepr(data_type="molecule", model_name="unimolv1", remove_hs=False,
                       pretrained_model_path=str(args.unimol_weights.resolve()),
                       pretrained_dict_path=str(args.unimol_dict.resolve()),
                       batch_size=8, use_cuda=args.device == "cuda")
    model.model.eval().requires_grad_(False)
    model.params.update(multi_process=False, seed=42)
    matrix = np.full((len(parents), WIDTHS["UniMol"]), np.nan, dtype=np.float32)
    failures = {}

    @torch.inference_mode()
    def batch(indices):
        smiles = parents.iloc[indices].drug_parent_smiles.tolist()
        # The official API returns only molecular CLS vectors in this mode.
        values = model.get_repr(smiles, return_atomic_reprs=False)
        matrix[indices] = checked_vectors("UniMol", np.asarray(values, dtype=np.float32), len(indices))

    for start in range(0, len(parents), 16):
        indices = list(range(start, min(start + 16, len(parents))))
        try:
            batch(indices)
        except (ValueError, RuntimeError, TypeError, IndexError):
            for index in indices:
                try:
                    batch([index])
                except (ValueError, RuntimeError, TypeError, IndexError) as error:
                    failures[parents.iloc[index].drug_parent_inchikey] = str(error)
    np.save(directory / "UniMol.npy", matrix)
    return dict(unimol_tools=version("unimol_tools"), torch=torch.__version__,
                rdkit=rdBase.rdkitVersion, model="unimolv1", readout="cls_repr",
                device=str(model.device), remove_hs=False, batch_size=8, seed=42, multi_process=False, failures=failures,
                checkpoint=dict(file=args.unimol_weights.name, sha256=sha256(args.unimol_weights)),
                dictionary=dict(file=args.unimol_dict.name, sha256=sha256(args.unimol_dict)))


def freeze(args):
    sys.path.insert(0, str(args.dose_api.resolve()))
    from tsg_round1.inputs import Inputs

    inputs = Inputs(args.inputs)
    coverage = fold_coverage(inputs)
    if inputs.source != "tenv" or sha256(args.kpgt_weights) != KPGT_BASE_SHA256:
        raise ValueError("requires accepted tenv inputs and the frozen KPGT base checkpoint")
    directory = args.output.resolve()
    directory.mkdir(parents=True)
    parents = inputs.parents[PARENT_COLUMNS].reset_index(drop=True)
    parents.to_csv(directory / "parents.csv", index=False)
    provenance = circular(parents, directory)
    if not np.array_equal(np.load(directory / "ECFP4C.npy"), inputs.morgan):
        raise ValueError("ECFP4C differs from the historical ECFP4 contract")
    subprocess.run([str(args.kpgt_python), str(Path(__file__).with_name("extract_kpgt.py")),
                    "--source", str(args.kpgt_source), "--weights", str(args.kpgt_weights),
                    "--parents", str(directory / "parents.csv"), "--output", str(directory),
                    "--device", args.device], check=True)
    kpgt_metadata = directory / "kpgt.json"
    provenance["KPGT"] = json.loads(kpgt_metadata.read_text())
    kpgt_metadata.unlink()
    provenance["KPGT"]["checkpoint_sha256"] = KPGT_BASE_SHA256
    kpgt_source_files = (
        "src/data/collator.py", "src/data/featurizer.py", "src/model/light.py", "src/model_config.py",
    )
    provenance["KPGT"]["source_files"] = {
        name: sha256(args.kpgt_source / name) for name in kpgt_source_files
    }
    provenance["E3FP"] = e3fp(parents, directory)
    provenance["UniMol"] = unimol(parents, directory, args)
    gaps = {}
    for arm in ARMS:
        matrix = np.load(directory / f"{arm}.npy", mmap_mode="r")
        missing = set(parents.loc[~np.isfinite(matrix).all(axis=1), "drug_parent_inchikey"])
        if missing:
            gaps[arm] = [
                dict(fold=fold, partition=part, missing_parents=sorted(
                    missing & set(inputs.obs.iloc[inputs.rows("cold-drug", fold, part)].drug_parent_inchikey)))
                for fold in range(3) for part in ("train", "validation", "test")
            ]
        else:
            checked_vectors(arm, matrix, len(parents))
    if gaps:
        write_json(directory / "coverage_gaps.json", dict(gaps=gaps, provenance=provenance))
        raise ValueError(f"representation coverage incomplete; see {directory / 'coverage_gaps.json'}")
    manifest = dict(schema=SCHEMA, widths=WIDTHS, coverage=coverage,
                    inputs_sha256=sha256(inputs.directory / "manifest.json"),
                    provenance=provenance, dtype="float32", feature_scaling="native",
                    files={path.name: sha256(path) for path in sorted(directory.iterdir())},
                    extraction_sha256=sha256(__file__),
                    kpgt_extraction_sha256=sha256(Path(__file__).with_name("extract_kpgt.py")))
    write_json(directory / "manifest.json", manifest)
    RepresentationInputs(inputs, directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("inputs", "output", "dose-api", "kpgt-python", "kpgt-source", "kpgt-weights",
                  "unimol-weights", "unimol-dict"):
        parser.add_argument(f"--{field}", required=True, type=Path)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    freeze(parser.parse_args())
