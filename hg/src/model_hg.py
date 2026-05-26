import torch
from torch_geometric.nn import SAGEConv, HeteroConv
import numpy as np

from src.utils_hg import get_batch_drug_subgraph


class HeteroGNN(torch.nn.Module):
    def __init__(self, node_input_dim=512, hidden_dim=128, out_dim=128, num_layers=2, include_edges=None):
        super().__init__()

        self.transforme_gene_node = torch.nn.Linear(node_input_dim, hidden_dim)
        self.transforme_drug_node = torch.nn.Linear(node_input_dim, hidden_dim)

        all_conv_layers = {
            ('gene', 'PPI', 'gene'): SAGEConv((-1, -1), hidden_dim),
            ('gene', 'GSEA', 'gene'): SAGEConv((-1, -1), hidden_dim),
            ('gene', 'PCC', 'gene'): SAGEConv((-1, -1), hidden_dim),
            ('drug', 'DDS', 'drug'): SAGEConv((-1, -1), hidden_dim),
            ('drug', 'DTI', 'gene'): SAGEConv((-1, -1), hidden_dim),
            ('gene', 'DTI', 'drug'): SAGEConv((-1, -1), hidden_dim)
        }

        if include_edges is not None:
            conv_layers = {}
            for edge_type, conv in all_conv_layers.items():
                if edge_type[1] in include_edges:
                    conv_layers[edge_type] = conv
                elif edge_type[1] == 'DTI' and 'DTI' in include_edges:
                    conv_layers[('gene', 'DTI', 'drug')] = SAGEConv((-1, -1), hidden_dim)
        else:
            conv_layers = all_conv_layers

        self.convs = torch.nn.ModuleList([HeteroConv(conv_layers) for _ in range(num_layers)])

        self.relation_predictors = torch.nn.ModuleDict({
            'DTI': torch.nn.Linear(hidden_dim * 2, 1),
            'DDS': torch.nn.Linear(hidden_dim * 2, 1),
            'PPI': torch.nn.Linear(hidden_dim * 2, 1)
        })

    def forward(self, x_dict, edge_index_dict):
        x_dict['gene'] = self.transforme_gene_node(x_dict['gene'])
        x_dict['drug'] = self.transforme_drug_node(x_dict['drug'])
        for conv in self.convs:
            x_dict = {key: x.relu() for key, x in conv(x_dict, edge_index_dict).items()}
        return x_dict

    def predict_bipartite_relations(self, src_emb, dst_emb, relation_type):
        num_src = src_emb.size(0)
        num_dst = dst_emb.size(0)

        src_expanded = src_emb.unsqueeze(1).expand(-1, num_dst, -1)
        dst_expanded = dst_emb.unsqueeze(0).expand(num_src, -1, -1)

        pair_repr = torch.cat([src_expanded, dst_expanded], dim=-1)
        pair_repr = pair_repr.view(-1, pair_repr.size(-1))

        scores = self.relation_predictors[relation_type](pair_repr)
        scores = scores.view(num_src, num_dst)

        return scores

    def predict_pairwise_relations(self, node_emb, relation_type):
        num_nodes = node_emb.size(0)

        src_expanded = node_emb.unsqueeze(1).expand(-1, num_nodes, -1)
        dst_expanded = node_emb.unsqueeze(0).expand(num_nodes, -1, -1)

        pair_repr = torch.cat([src_expanded, dst_expanded], dim=-1)
        pair_repr = pair_repr.view(-1, pair_repr.size(-1))

        scores = self.relation_predictors[relation_type](pair_repr)
        scores = scores.view(num_nodes, num_nodes)

        return scores

    def predict_relations(self, x_dict, edge_type=None):
        drug_emb = x_dict['drug']
        gene_emb = x_dict['gene']

        if edge_type == 'DTI':
            scores = self.predict_bipartite_relations(drug_emb, gene_emb, 'DTI')
        elif edge_type == 'DDS':
            scores = self.predict_pairwise_relations(drug_emb, 'DDS')
        elif edge_type == 'PPI':
            scores = self.predict_pairwise_relations(gene_emb, 'PPI')
        else:
            raise ValueError(f"Invalid edge type: {edge_type}")

        return scores


def compute_relation_matrix(drug_embedding, gene_embedding, similarity_metric='cosine'):
    if similarity_metric == 'cosine':
        norm_drug = np.linalg.norm(drug_embedding, axis=1, keepdims=True)
        norm_gene = np.linalg.norm(gene_embedding, axis=1, keepdims=True)
        drug_normalized = np.divide(drug_embedding, norm_drug, where=norm_drug!=0)
        gene_normalized = np.divide(gene_embedding, norm_gene, where=norm_gene!=0)
        relation_matrix = np.dot(drug_normalized, gene_normalized.T)
    elif similarity_metric == 'dot':
        relation_matrix = np.dot(drug_embedding, gene_embedding.T)
    else:
        raise ValueError("Similarity metric must be 'cosine' or 'dot'")

    return relation_matrix


def predict_all_dti_scores(model, data, device, logger, batch_size=128):
    model.eval()

    num_drugs = data.x_dict['drug'].shape[0]
    num_genes = data.x_dict['gene'].shape[0]

    scores_matrix = np.zeros((num_drugs, num_genes))

    with torch.no_grad():
        for i in range(0, num_drugs, batch_size):
            batch_end = min(i + batch_size, num_drugs)
            drug_indices = torch.arange(i, batch_end).to(device)

            logger.info(f'Processing drugs batch from index {i} to {batch_end-1}')

            subgraph = get_batch_drug_subgraph(data, drug_indices, device)
            subgraph = subgraph.to(device)

            x_dict = model(subgraph.x_dict, subgraph.edge_index_dict)

            batch_scores = model.predict_relations(x_dict, edge_type='DTI')

            batch_scores = batch_scores.cpu().numpy()
            batch_scores = 1 / (1 + np.exp(-batch_scores))

            scores_matrix[i:batch_end] = batch_scores

            del subgraph, x_dict, batch_scores
            torch.cuda.empty_cache()

    return scores_matrix
