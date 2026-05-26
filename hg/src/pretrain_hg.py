import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
import os.path as osp
import time
import logging
import argparse
import random

import torch
import numpy as np
from torch_geometric.data import HeteroData
from torch_geometric.loader import NeighborLoader

from src.utils_hg import log_nested_dict, generate_multi_relation_samples, \
                ContrastiveLoss, get_embeddings, get_batch_drug_subgraph
from src.model_hg import HeteroGNN, compute_relation_matrix, predict_all_dti_scores


HG_ROOT = osp.dirname(osp.dirname(osp.abspath(__file__)))


def get_output_dir(timestamp):
    return osp.join(HG_ROOT, 'model', timestamp)


def set_random_seed(seed=10):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)


class EarlyStopping():

    def __init__(self, mode='higher', patience=None, metric=None, n_fold=None, folder=None, logger=None):

        filepath = os.path.join(folder, '{}_fold_early_stop.pth'.format(n_fold))

        if metric is not None:
            assert metric in ['r2', 'mae', 'rmse', 'roc_auc_score', 'pr_auc_score', 'mse'], \
                "Expect metric to be 'r2' or 'mae' or " \
                "'rmse' or 'roc_auc_score' or 'mse', got {}".format(metric)
            if metric in ['r2', 'roc_auc_score', 'pr_auc_score']:
                mode = 'higher'
            if metric in ['mae', 'rmse', 'mse']:
                mode = 'lower'

        assert mode in ['higher', 'lower']
        self.mode = mode
        if self.mode == 'higher':
            self._check = self._check_higher
        else:
            self._check = self._check_lower

        self.patience = patience
        self.counter = 0
        self.filepath = filepath
        self.best_score = None
        self.early_stop = False
        self.logger = logger

    def _check_higher(self, score, prev_best_score):
        return score > prev_best_score

    def _check_lower(self, score, prev_best_score):
        return score < prev_best_score

    def step(self, score, model, current_epoch, optimizer):
        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(model, current_epoch, optimizer)
        elif self._check(score, self.best_score):
            self.best_score = score
            self.save_checkpoint(model, current_epoch, optimizer)
            self.counter = 0
        else:
            self.counter += 1
            self.logger.info(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True

        return self.early_stop

    def save_checkpoint(self, model, current_epoch, optimizer):
        torch.save({
            'epoch': current_epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
        }, self.filepath)

    def load_checkpoint(self, model, optimizer):
        checkpoint = torch.load(self.filepath)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        return model


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--node_input_dim', type=int, default=512)
    parser.add_argument('--hidden_dim', type=int, default=256)
    parser.add_argument('--out_dim', type=int, default=256)
    parser.add_argument('--num_layers', type=int, default=3)
    parser.add_argument('--learning_rate', type=float, default=0.0005)
    parser.add_argument('--hg_path', type=str, default=None)
    parser.add_argument('--output_dir', type=str, default=None)
    parser.add_argument('--include_edges', type=str, nargs='+', default=['PPI', 'DTI', 'DDS'])
    parser.add_argument('--seed', type=int, default=4242)
    parser.add_argument('--device', type=str, default='cuda:0')

    parser.add_argument('--epochs', type=int, default=300)
    parser.add_argument('--patience', type=int, default=30)
    parser.add_argument('--num_neg_samples', type=int, default=5)
    parser.add_argument('--num_neighbors', nargs='+', type=int, default=[35, 20, 10])
    parser.add_argument('--temperature', type=float, default=0.5)
    parser.add_argument('--batch_size', type=int, default=4)
    parser.add_argument('--global_saved_model', type=str, default=None)
    parser.add_argument('--num_iters', type=int, default=2)
    parser.add_argument('--output_relation_matrix', type=bool, default=False)
    parser.add_argument('--infer_dti_scores', type=bool, default=False)
    parser.add_argument('--conseq_patience', type=int, default=3)
    parser.add_argument('--turn_patience', type=int, default=10)

    return parser.parse_args()



def train_contrastive_batch(model, data, optimizer, loader, epochs=100, patience=30, num_neg_samples=5, temperature=0.5, expt_folder=None, logger=None, nfold='global', conseq_patience=3, turn_patience=10):

    model.train()
    contrastive_loss_fn = ContrastiveLoss(temperature)
    stopper = EarlyStopping(mode='lower', metric='mse', patience=patience, n_fold=nfold, folder=expt_folder, logger=logger)

    prev_loss = float('inf')
    conseq_increase = 0
    total_increase = 0

    for epoch in range(epochs):
        epoch_loss = 0

        i = 0
        for batch in loader:
            optimizer.zero_grad()
            i+=1

            x_dict = model(batch.x_dict, batch.edge_index_dict)

            contrastive_loss = torch.tensor(0.0).to(batch.x_dict['drug'].device)
            pos_samples, neg_samples = generate_multi_relation_samples(batch, num_neg_samples)
            if not pos_samples or not neg_samples:
                logger.warning("Empty samples generated, skipping batch")
                continue

            valid_edges = 0
            for edge_type in pos_samples:
                try:
                    pos_src_idx = pos_samples[edge_type][0]
                    pos_dst_idx = pos_samples[edge_type][1]
                    neg_dst_idx = neg_samples[edge_type][1]

                    z_i = x_dict[edge_type[0]][pos_src_idx]
                    z_j = x_dict[edge_type[2]][pos_dst_idx]
                    z_neg = x_dict[edge_type[2]][neg_dst_idx]

                    edge_loss = contrastive_loss_fn(z_i, z_j, z_neg)
                    if not torch.isnan(edge_loss):
                        contrastive_loss += edge_loss
                        valid_edges += 1

                except Exception as e:
                    logger.error(f"Error processing edge type {edge_type}: {str(e)}")
                    continue

            logger.info(f'batch_contrast_loss:{contrastive_loss}')

            loss = contrastive_loss
            logger.info(f'batch_loss:{loss}')

            if loss != 0.:
                loss.backward()
                optimizer.step()

            epoch_loss += loss.item()

        avg_loss = epoch_loss / len(loader)

        if avg_loss > prev_loss:
            conseq_increase += 1
            total_increase += 1
        else:
            conseq_increase = 0
        prev_loss = avg_loss

        logger.info(f'Epoch {epoch:03d}, '
                   f'Loss: {avg_loss:.4f}, '
                   f'conseq_increase: {conseq_increase}, '
                   f'total_increase: {total_increase}')

        if conseq_increase >= conseq_patience:
            logger.info(f'Consecutive early stopping at epoch {epoch}: loss increased for {conseq_patience} consecutive epochs')
            break

        if total_increase >= turn_patience:
            logger.info(f'Turn early stopping at epoch {epoch}: loss increased for the {turn_patience}th time')
            break

        early_stop = stopper.step(avg_loss, model, epoch, optimizer)
        if early_stop:
            logger.info(f'Early stopping at epoch {epoch}')
            break

    best_model = stopper.load_checkpoint(model, optimizer)
    return best_model




def main():

    start_time = time.time()
    timestamp = time.strftime('%Y%m%d_%H%M%S', time.localtime())

    args = get_args()

    if args.hg_path is None:
        args.hg_path = osp.join(HG_ROOT, 'data') + '/'

    if args.output_dir:
        expt_folder = osp.join(args.output_dir, 'logs', timestamp)
    else:
        expt_folder = str(get_output_dir(timestamp))

    if not os.path.exists(expt_folder):
        os.makedirs(expt_folder)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(expt_folder+f'/{timestamp}.log')
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    logger.info('\n---------args-----------')
    log_nested_dict(vars(args), logger)
    logger.info('\n')

    set_random_seed(args.seed)
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    model = HeteroGNN(
        hidden_dim=args.hidden_dim,
        out_dim=args.out_dim,
        num_layers=args.num_layers,
        include_edges=args.include_edges
    ).to(device)
    logger.info(model)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)

    logger.info("Loading data...")
    gene_feat = torch.tensor(np.load(args.hg_path + 'all_gene_node_feat_19392.npy'), dtype=torch.float)
    drug_feat = torch.tensor(np.load(args.hg_path + 'all_drug_node_feat_8981.npy'), dtype=torch.float)

    data = HeteroData()
    data['gene'].x = gene_feat
    data['drug'].x = drug_feat

    edge_types = ['PPI', 'GSEA', 'PCC', 'DTI', 'DDS']
    edge_files = [
        'PPI_edge_index_all_901260_pairs_weighted.npy',
        'GSEA_edge_indexes_all_pairs_426904_weighted.npy',
        'PCC_edge_Thr04_18674_pairs.npy',
        'DTI_edge_index_all_12890_pairs.npy',
        'DDS_edge_index_all_287834_pairs_weighted.npy'
    ]

    for edge_type, edge_file in zip(edge_types, edge_files):
        if edge_type in args.include_edges:
            edge_data = np.load(args.hg_path + edge_file, allow_pickle=True).item()
            edge_index = torch.tensor(edge_data['edge_index'], dtype=torch.long)
            edge_attr = torch.tensor(edge_data['edge_attr'], dtype=torch.float) if 'weighted' in edge_file else None

            if edge_type in ['PPI', 'GSEA']:
                data['gene', edge_type, 'gene'].edge_index = edge_index
                data['gene', edge_type, 'gene'].edge_attr = edge_attr
            elif edge_type == 'PCC':
                data['gene', edge_type, 'gene'].edge_index = edge_index
            elif edge_type == 'DTI':
                data['drug', edge_type, 'gene'].edge_index = edge_index
                data['gene', edge_type, 'drug'].edge_index = torch.stack([edge_index[1], edge_index[0]], dim=0)
            elif edge_type == 'DDS':
                data['drug', edge_type, 'drug'].edge_index = edge_index
                data['drug', edge_type, 'drug'].edge_attr = edge_attr
            else:
                logger.info(f'{edge_type} edge_type is not assigned!')


    logger.info('Checking isolated nodes...')
    all_drug_indices = set(range(data['drug'].x.size(0)))
    all_gene_indices = set(range(data['gene'].x.size(0)))
    connected_drugs = set()
    connected_genes = set()

    for edge_type, edge_index_dict in data.edge_index_dict.items():
        src_type, relation, dst_type = edge_type
        edge_index = edge_index_dict
        if src_type == 'drug':
            connected_drugs.update(edge_index[0].cpu().numpy())
        elif src_type == 'gene':
            connected_genes.update(edge_index[0].cpu().numpy())

        if dst_type == 'drug':
            connected_drugs.update(edge_index[1].cpu().numpy())
        elif dst_type == 'gene':
            connected_genes.update(edge_index[1].cpu().numpy())

    isolated_drugs = all_drug_indices - connected_drugs
    isolated_genes = all_gene_indices - connected_genes

    logger.info(f"Total number of drug nodes: {len(all_drug_indices)}")
    logger.info(f"Number of isolated drug nodes: {len(isolated_drugs)}")
    logger.info(f"Total number of gene nodes: {len(all_gene_indices)}")
    logger.info(f"Number of isolated gene nodes: {len(isolated_genes)}")

    if len(isolated_drugs) > 0 or len(isolated_genes) > 0:
        np.save(f"{expt_folder}/isolated_drug_indices.npy", np.array(list(isolated_drugs)))
        np.save(f"{expt_folder}/isolated_gene_indices.npy", np.array(list(isolated_genes)))

    data = data.to(device)

    if args.global_saved_model:
        logger.info(f"Loading pretrained global model from {args.global_saved_model}")
        model.load_state_dict(torch.load(args.global_saved_model)['model_state_dict'])

        logger.info('loading Phase 1 embeddings...')
        drug_embeddings, gene_embeddings = get_embeddings(model, data)
        phase1_drug_embeddings = drug_embeddings.cpu().detach().numpy()
        phase1_gene_embeddings = gene_embeddings.cpu().detach().numpy()
    else:
        logger.info("Starting Phase 1: Global Graph Training...")
        loader = NeighborLoader(
            data,
            num_neighbors=args.num_neighbors,
            batch_size=args.batch_size,
            input_nodes=('drug', torch.arange(data['drug'].x.size(0)))
        )

        model = train_contrastive_batch(
            model=model,
            data=data,
            optimizer=optimizer,
            loader=loader,
            epochs=args.epochs,
            patience=args.patience,
            num_neg_samples=args.num_neg_samples,
            temperature=args.temperature,
            expt_folder=expt_folder,
            logger=logger,
            nfold='global',
            conseq_patience=args.conseq_patience,
            turn_patience=args.turn_patience
        )

        logger.info(f'Phase 1 completed! Training time: {int(time.time() - start_time)/60:.2f}min')

        num_dti_edges = sum(len(edge_index[0]) for edge_type, edge_index in data.edge_index_dict.items() if edge_type[1] == 'DTI')
        logger.info(f'Number of DTI edges in data: {num_dti_edges}')

        logger.info('Saving Phase 1 embeddings...')
        drug_embeddings, gene_embeddings = get_embeddings(model, data)
        phase1_drug_embeddings = drug_embeddings.cpu().detach().numpy()
        phase1_gene_embeddings = gene_embeddings.cpu().detach().numpy()
        np.save(f"{expt_folder}/phase1_drug_embeddings.npy", phase1_drug_embeddings)
        np.save(f"{expt_folder}/phase1_gene_embeddings.npy", phase1_gene_embeddings)

        embeddings_dir = osp.join(HG_ROOT, 'embeddings')
        np.save(osp.join(embeddings_dir, "HG_drug_embeddings.npy"), phase1_drug_embeddings)
        logger.info(f"HG drug embeddings also saved to {osp.join(embeddings_dir, 'HG_drug_embeddings.npy')}")

    if args.output_relation_matrix:
        logger.info('Extracting relation matrix...')
        phase1_gene_embeddings_mean = np.mean(phase1_gene_embeddings, axis=0, keepdims=True)
        phase1_gene_embeddings_arr = np.concatenate([phase1_gene_embeddings, phase1_gene_embeddings_mean], axis=0)
        logger.info('loadding gene_978_idx...')
        gene_978_idx = np.load(args.hg_path + '978_gene_node_idx.npy', allow_pickle=True)
        phase1_gene_embeddings_arr_978 = phase1_gene_embeddings_arr[gene_978_idx ,:]
        relation_matrix = compute_relation_matrix(phase1_drug_embeddings, phase1_gene_embeddings_arr_978, similarity_metric='cosine')
        np.save(f"{expt_folder}/phase1_relation_matrix.npy", relation_matrix)

    if args.infer_dti_scores:
        logger.info('Inferring DTI scores...')
        dti_scores_matrix = predict_all_dti_scores(model, data, device, logger, batch_size=64)
        logger.info(f'dti_scores_matrix.shape:{dti_scores_matrix.shape}')
        phase1_dti_scores_mean = np.mean(dti_scores_matrix, axis=-1, keepdims=True)
        phase1_dti_arr = np.concatenate([dti_scores_matrix, phase1_dti_scores_mean], axis=-1)
        logger.info('loadding gene_978_idx...')
        gene_978_idx = np.load(args.hg_path + '978_gene_node_idx.npy', allow_pickle=True)
        phase1_dti_arr_978 = phase1_dti_arr[:, gene_978_idx]
        logger.info(f'phase1_dti_arr_978.shape:{phase1_dti_arr_978.shape}')
        np.save(f"{expt_folder}/phase1_dti_scores_matrix_978.npy", phase1_dti_arr_978)


    logger.info(f'Total training time: {int(time.time() - start_time)/60:.2f}min')
    logger.info("Training completed successfully!")

if __name__ == '__main__':
    main()
