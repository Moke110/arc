"""Use the frozen KPGT source and checkpoint without writing to source/."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def extract(args):
    sys.path.insert(0, str(args.source.resolve()))
    import dgl
    import scipy
    from rdkit import Chem, rdBase
    import torch
    from src.data.collator import Collator_tune
    from src.data.descriptors.rdNormalizedDescriptors import RDKit2DNormalized
    from src.data.featurizer import N_ATOM_TYPES, N_BOND_TYPES, Vocab, smiles_to_graph_tune
    from src.model.light import LiGhTPredictor
    from src.model_config import config_dict

    torch.set_num_threads(2)
    config = config_dict["base"]
    options = {name: config[name] for name in (
        "d_node_feats", "d_edge_feats", "d_g_feats", "d_hpath_ratio", "n_mol_layers",
        "path_length", "n_heads", "n_ffn_dense_layers",
    )}
    model = LiGhTPredictor(**options, n_node_types=Vocab(N_ATOM_TYPES, N_BOND_TYPES).vocab_size,
                          input_drop=0, attn_drop=0, feat_drop=0)
    weights = torch.load(args.weights, map_location="cpu", weights_only=True)
    model.load_state_dict({k[7:] if k.startswith("module.") else k: value for k, value in weights.items()})
    model.to(args.device).eval().requires_grad_(False)
    descriptors = RDKit2DNormalized()
    collate = Collator_tune(config["path_length"])
    parents = pd.read_csv(args.parents)
    results = np.full((len(parents), 2304), np.nan, dtype=np.float32)
    failures = {}

    def data(smiles):
        graph = smiles_to_graph_tune(smiles, max_length=config["path_length"], n_virtual_nodes=2)
        if graph is None:
            raise ValueError("KPGT graph conversion failed")
        molecule = Chem.MolFromSmiles(smiles)
        fp = list(Chem.RDKFingerprint(molecule, minPath=1, maxPath=7, fpSize=512))
        md = descriptors.process(smiles)[1:]
        return (smiles, graph, torch.tensor(fp, dtype=torch.float32),
                torch.tensor(np.where(np.isnan(md), 0, md), dtype=torch.float32), torch.zeros(1))

    @torch.inference_mode()
    def batch(indices):
        samples = [data(parents.iloc[i].drug_parent_smiles) for i in indices]
        _, graph, fp, md, _ = collate(samples)
        values = model.generate_fps(graph.to(args.device), fp.to(args.device), md.to(args.device))
        values = values.cpu().numpy()
        if values.shape != (len(indices), 2304) or not np.isfinite(values).all():
            raise ValueError("native KPGT readout must return 2304 finite values per parent")
        results[indices] = values

    smoke = list(range(min(3, len(parents))))
    batch(smoke)
    print(f"KPGT native smoke passed for {len(smoke)} parents, shape {results[smoke].shape}",
          flush=True)
    for start in range(len(smoke), len(parents), 8):
        indices = list(range(start, min(start + 8, len(parents))))
        try:
            batch(indices)
        except (ValueError, RuntimeError, TypeError) as error:
            print(f"KPGT batch failed, retrying individually: {error}", file=sys.stderr)
            for i in indices:
                try:
                    batch([i])
                except (ValueError, RuntimeError, TypeError) as error:
                    failures[parents.iloc[i].drug_parent_inchikey] = str(error)

    np.save(args.output / "KPGT.npy", results)
    info = dict(device=args.device, torch=torch.__version__, dgl=dgl.__version__,
                rdkit=rdBase.rdkitVersion, numpy=np.__version__, scipy=scipy.__version__,
                smoke_parents=parents.iloc[smoke].drug_parent_inchikey.tolist(),
                smoke_shape=list(results[smoke].shape), failures=failures, config=options,
                readout="LiGhTPredictor.generate_fps", batch_size=8,
                descriptors="RDKit2DNormalized", fingerprint=dict(minPath=1, maxPath=7, fpSize=512),
                virtual_nodes=2)
    (args.output / "kpgt.json").write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("source", "weights", "parents", "output"):
        parser.add_argument(f"--{field}", required=True, type=Path)
    parser.add_argument("--device", choices=("cpu", "cuda"), required=True)
    extract(parser.parse_args())
