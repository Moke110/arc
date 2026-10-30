#!/usr/bin/env python
"""CMonge setup smoke: official demo pipeline with local paths and a short run."""

import sys

# Keep `pert_model/*/source/` byte-identical to its download (ADR 0030): this
# smoke imports model code from the vendored tree, so byte-code must not land
# there. Set before any other import; the run.sh wrapper exports the same flag.
sys.dont_write_bytecode = True

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def find_repo_root(start: Path) -> Path:
    """Walk up to the directory containing .git (ADR 0030 §4)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"cannot find the repository root above {start}")


ROOT = find_repo_root(HERE)
MODEL = ROOT / "pert_model" / "CMonge" / "source"
OUT = HERE / "output"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(MODEL))

from cmonge.datasets.conditional_loader import ConditionalDataModule  # noqa: E402
from cmonge.trainers.ae_trainer import AETrainerModule  # noqa: E402
from cmonge.trainers.conditional_monge_trainer import ConditionalMongeTrainer  # noqa: E402
from cmonge.utils import load_config  # noqa: E402


def main():
    config = load_config(str(MODEL / "configs" / "demo_config.yml"))
    config.data.file_path = str(MODEL / "data" / "dummy_data.h5ad")
    config.data.ae_config_path = str(MODEL / "configs" / "autoencoder-demo.yml")
    config.model.embedding.smile_path = str(MODEL / "data" / "lincs_trapnell.smiles")
    config.model.embedding.drug_to_smile_path = str(MODEL / "data" / "trapnell_drugs_smiles.csv")
    config.model.embedding.model_dir = str(OUT / "embed")
    config.model.num_train_iters = 10
    config.ae.training.model_dir = str(OUT / "ae")
    config.ae.training.n_epochs = 2

    config.data.ae = True
    config.data.reduction = None
    datamodule = ConditionalDataModule(config.data, config.condition)
    ae_trainer = AETrainerModule(config.ae)
    ae_trainer.train(datamodule)

    config.data.ae = False
    config.ae.model.act_fn = "gelu"
    config.data.reduction = "ae"
    datamodule = ConditionalDataModule(
        config.data, config.condition, ae_config=config.ae
    )
    trainer = ConditionalMongeTrainer(
        jobid=1,
        logger_path=str(OUT / "demo_logs.yml"),
        config=config.model,
        datamodule=datamodule,
    )
    trainer.train(datamodule)
    trainer.evaluate(datamodule)

    meta = {
        "model": "CMonge",
        "n_genes": 100,
        "conditions": list(config.condition.conditions),
        "num_train_iters": config.model.num_train_iters,
        "ae_epochs": config.ae.training.n_epochs,
        "status": "ok",
        "note": "official demo data; dose is a categorical condition label",
    }
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
