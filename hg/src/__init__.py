def __getattr__(name):
    if name == 'HeteroGNN':
        from src.model_hg import HeteroGNN
        return HeteroGNN
    elif name in ('ContrastiveLoss', 'DrugSpecificLoss', 'get_embeddings', 'generate_multi_relation_samples', 'get_batch_drug_subgraph', 'get_single_drug_subgraph', 'update_graph_with_predictions', 'get_high_confidence_edges'):
        from src.utils_hg import ContrastiveLoss, DrugSpecificLoss, get_embeddings, generate_multi_relation_samples, get_batch_drug_subgraph, get_single_drug_subgraph, update_graph_with_predictions, get_high_confidence_edges
        return locals()[name]
    elif name in ('compute_relation_matrix', 'predict_all_dti_scores'):
        from src.model_hg import compute_relation_matrix, predict_all_dti_scores
        return compute_relation_matrix if name == 'compute_relation_matrix' else predict_all_dti_scores
    raise AttributeError(f"module 'src' has no attribute {name}")
