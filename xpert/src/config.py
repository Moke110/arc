import os
import yaml
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = PROJECT_ROOT.parent


def load_config(config_name: str) -> dict:
    config_path = PROJECT_ROOT / f"{config_name}.yaml"
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    config = resolve_paths(config)
    return config


def resolve_paths(config: dict) -> dict:
    xpert_root = str(PROJECT_ROOT)
    workspace_root = str(WORKSPACE_ROOT)

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
