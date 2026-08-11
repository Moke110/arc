import os
import yaml
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = PROJECT_ROOT.parent


def _resolve_config_path(config_name: str) -> Path:
    """Resolve a config name to an existing YAML file path.

    Order of resolution:
      1. absolute path as given (with .yaml appended if missing)
      2. path relative to the current working directory
      3. path relative to the XPert project root (backward-compatible)
    """
    candidate = Path(config_name)
    if candidate.suffix != ".yaml":
        candidate = candidate.with_suffix(".yaml")

    if candidate.is_absolute():
        return candidate

    cwd_candidate = Path.cwd() / candidate
    if cwd_candidate.exists():
        return cwd_candidate

    return PROJECT_ROOT / candidate


def get_config_root(config_name: str) -> Path:
    """Return the project root that owns the given config file.

    Configs live at <root>/<stage>/config/<name>.yaml, so the root is three
    directories up from the config file (config -> <stage> -> <root>). This lets
    a workflow own its configs (e.g. workflow/xpert-fda-miri-v1/05_xpert_finetune/)
    while the model library resolves paths against the owning root.
    """
    config_path = _resolve_config_path(config_name)
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_name}")
    return config_path.resolve().parent.parent.parent


def load_config(config_name: str) -> dict:
    config_path = _resolve_config_path(config_name)
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    config_root = config_path.resolve().parent.parent.parent
    config = resolve_paths(config, config_root)
    return config


def resolve_paths(config: dict, config_root: Path = PROJECT_ROOT) -> dict:
    xpert_root = str(config_root)

    dataset = config.get('dataset', {})
    for key in ['l1000_sdst_data_root', 'l1000_mdmt_data_root',
                'l1000_mdmt_pretrain_data_root', 'l1000_mdmt_full_data_root',
                'l1000_mdmt_example_data_root',
                'CDSDB_data_root', 'PANACEA_data_root', 'DTOX_data_root']:
        if key in dataset and dataset[key] is not None:
            if not os.path.isabs(dataset[key]):
                dataset[key] = os.path.join(xpert_root, dataset[key])

    model = config.get('model', {})
    attn = model.get('ATTN', {})
    if 'ppi_gene_vector_path' in attn and attn['ppi_gene_vector_path'] is not None:
        if not os.path.isabs(attn['ppi_gene_vector_path']):
            attn['ppi_gene_vector_path'] = os.path.join(xpert_root, attn['ppi_gene_vector_path'])

    drug_embeddings = config.get('drug_embeddings', {})
    for key, path in drug_embeddings.items():
        if path is not None and not os.path.isabs(path):
            drug_embeddings[key] = os.path.join(xpert_root, path)

    return config


def get_output_dir(stage: str, timestamp: str, dataset: str = '', split_method: str = '') -> Path:
    if stage == 'pretrain':
        return PROJECT_ROOT / 'pretrain' / 'output' / f'{dataset}_{split_method}' / timestamp
    elif stage == 'finetune':
        return PROJECT_ROOT / 'finetune' / 'output' / f'{dataset}_{split_method}' / timestamp
    elif stage == 'infer':
        return PROJECT_ROOT / 'metrics' / 'output' / timestamp
    else:
        return PROJECT_ROOT / 'experiment' / timestamp
