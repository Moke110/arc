import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

import json
import torch
from torch.amp import GradScaler, autocast
import numpy as np
import pandas as pd
import argparse
import time
import os.path as osp
import logging
from tqdm import tqdm
import pickle

from src.config import load_config, get_config_root, get_output_dir, PROJECT_ROOT
from src.utils import EarlyStopping, save_loss_fig, set_random_seed, log_nested_dict
from src.utils import load_dataloader, load_test_dataloader, load_infer_dataloader
from src.utils import mse_loss_ls_sum, pcc_loss_sum, cos_loss_sum, set_requires_grad
from src.metrics import get_metrics_new, get_metrics_infer
from src.model_xpert import XPertNet


def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', default='train', help='train or valid or test or infer')
    parser.add_argument('--split_method', default='random', help='Split method key from config split_methods. Pretrain: pretrain_split. Finetune: random, warm-drug, cold-drug.')
    parser.add_argument('--drug_feat', default='unimol', help='primary molecular feature: unimol, KPGT, morgan, smi')
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--model', default='XPert', help='XPert')
    parser.add_argument('--config', default='../../workflow/xpert-fda-miri-v1/03_xpert_pretrain/config/l1000', help='config path relative to project root (without .yaml)')
    parser.add_argument('--seed', default=2024, type=int)
    parser.add_argument('--dataset', default='l1000_sdst', help='l1000_sdst, l1000_mdmt, l1000_mdmt_pretrain, l1000_mdmt_full, cdsdb, panacea, dtoxs')
    parser.add_argument('--pretrain_dataset', type=str, default='', help='pretrain dataset name for finetune output dir naming')
    parser.add_argument('--pretrained_mode', default='global', help='global, specific')

    parser.add_argument('--include_cell_idx', type=bool, default=False)
    parser.add_argument('--wo_HG', type=bool, default=False)
    parser.add_argument('--wo_atom', type=bool, default=False)
    parser.add_argument('--wo_atom_HG', type=bool, default=False)
    parser.add_argument('--wo_unimol', type=bool, default=False)
    parser.add_argument('--wo_ppi', type=bool, default=False)
    parser.add_argument('--use_gene_pos_emed', type=bool, default=False)

    parser.add_argument('--use_gradscaler', type=bool, default=False)
    parser.add_argument('--lr_scheduler', type=bool, default=True)
    parser.add_argument('--resume_from', type=str, default=None)
    parser.add_argument('--saved_model_path', type=str, default=None)
    parser.add_argument('--saved_model', type=str, default=None)
    parser.add_argument('--pretrained_model', type=str, default=None)
    parser.add_argument('--output_profile', type=bool, default=False)
    parser.add_argument('--output_attention', type=bool, default=False)
    parser.add_argument('--output_cls_embed', type=bool, default=False)

    parser.add_argument('--weighted_loss', type=bool, default=False)
    parser.add_argument('--kl_loss', type=bool, default=False)

    # training hyperparameters (replaces config['train'])
    parser.add_argument('--train_lr', type=float, default=0.004)
    parser.add_argument('--weight_decay', type=float, default=1e-5)
    parser.add_argument('--batch_size', type=int, default=64)
    parser.add_argument('--num_epochs', type=int, default=2500)
    parser.add_argument('--init_epoch', type=int, default=70, help='warmup epochs before full loss backprop')
    parser.add_argument('--loss_weight', type=str, default='0.2,0.003,0.2,1', help='comma-separated: L_pert,L_ctl,L_deg,L_pcc')
    parser.add_argument('--patience', type=int, default=50)
    parser.add_argument('--bandwidth', type=float, default=0.3)

    # finetune hyperparameters (replaces config['finetune'])
    parser.add_argument('--finetune_lr', type=float, default=0.001)
    parser.add_argument('--finetune_weight_decay', type=float, default=1e-5)
    parser.add_argument('--finetune_params', type=str, default='all')
    parser.add_argument('--exclude_keys', type=str, default='')

    # drug embedding selection (resolved via config['drug_embeddings'])
    parser.add_argument('--drug_embeds', type=str, default='unimol,HG', help='comma-separated embedding names: unimol,KPGT,morgan,smi,HG')

    return parser.parse_args()



def train(model, opt, dataloader, args, config, scaler=None, epoch=0, loss_func=None, kl_loss=None):
    model.train()
    total_loss = 0.
    total_loss1 = 0.
    total_loss2 = 0.
    total_loss3 = 0.
    total_loss4 = 0.
    num_batches = 0

    mse_loss_func = mse_loss_ls_sum
    MultiClassLoss = torch.nn.CrossEntropyLoss(reduction='sum')
    a, b, c, d = args.loss_weight

    init_epoch = args.init_epoch

    for data in tqdm(dataloader, desc=f'Train epoch {epoch}', leave=False):

        model.zero_grad()

        if scaler:
            with autocast(device_type='cuda'):
                trt_output, ctl_output, deg_output, trt_raw_data, ctl_raw_data, _, cell_class_true, cell_class_predict = model(data)
                deg_true = trt_raw_data-ctl_raw_data
                num_samples =  trt_raw_data.shape[0]
                num_batches += num_samples
                loss1 = mse_loss_ls_sum(trt_output, trt_raw_data)
                if cell_class_predict is not None:
                    cell_class_predict_1, cell_class_predict_2 = cell_class_predict
                    loss2 = MultiClassLoss(cell_class_predict_1, cell_class_true) + MultiClassLoss(cell_class_predict_2, cell_class_true)
                else:
                    loss2 = mse_loss_ls_sum(ctl_output, ctl_raw_data) 
                
                loss3 = mse_loss_func(deg_output, deg_true) 
                loss4 = pcc_loss_sum(deg_output, deg_true)

                weighted_loss = loss1 * a + loss2 * b + loss3 * c + loss4 * d
                if cell_class_predict is not None:
                    batch_weighted_loss = torch.sqrt(loss1/num_samples) * a + (loss2/num_samples) * b + torch.sqrt(loss3/num_samples) * c + (loss4/num_samples) * d
                else:
                    batch_weighted_loss = torch.sqrt(loss1/num_samples) * a + torch.sqrt(loss2/num_samples) * b + torch.sqrt(loss3/num_samples) * c + (loss4/num_samples) * d
                total_loss += weighted_loss.item()
                total_loss1 += loss1.item()
                total_loss2 += loss2.item()
                total_loss3 += loss3.item()
                total_loss4 += loss4.item()

            if epoch < init_epoch:
                scaler.scale(batch_weighted_loss).backward()
            else:
                scaler.scale(weighted_loss).backward()
            
            scaler.step(opt)
            scaler.update()

        else:
            trt_output, ctl_output, deg_output, trt_raw_data, ctl_raw_data, _ = model(data)
            num_samples =  trt_raw_data.shape[0]
            num_batches += num_samples
            loss1 = mse_loss_ls_sum(trt_output, trt_raw_data)
            loss2 = mse_loss_ls_sum(ctl_output, ctl_raw_data)
            loss3 = mse_loss_ls_sum(deg_output, trt_raw_data-ctl_raw_data)
            loss4 = pcc_loss_sum(deg_output, trt_raw_data-ctl_raw_data)
            weighted_loss = loss1 * a + loss2 * b + loss3 * c + loss4 * d
            batch_weighted_loss = torch.sqrt(loss1/num_samples) * a + torch.sqrt(loss2/num_samples) * b + torch.sqrt(loss3/num_samples) * c + (loss4/num_samples) * d
            total_loss += weighted_loss.item()
            total_loss1 += loss1.item()
            total_loss2 += loss2.item()
            total_loss3 += loss3.item()
            total_loss4 += loss4.item()

            if epoch < init_epoch:
                batch_weighted_loss.backward()
            else:
                weighted_loss.backward()
            
            opt.step()


    avg_total_loss = total_loss / num_batches
    avg_loss1 = np.sqrt(total_loss1 / num_batches) if total_loss1 else -999
    if cell_class_predict is not None:
        avg_loss2 = (total_loss2 / num_batches) if total_loss2 else -999
    else:
        avg_loss2 = np.sqrt(total_loss2 / num_batches) if total_loss2 else -999
    avg_loss3 = np.sqrt(total_loss3 / num_batches) if total_loss3 else -999
    avg_loss4 = total_loss4 / num_batches if total_loss4 else -999

    return avg_total_loss, avg_loss1, avg_loss2, avg_loss3, avg_loss4




def validate(model, dataloader, args, config, flag=False, scaler=None, expt_folder=None, fold=None):
    model.eval()
    # Handle empty dataloader (e.g., dtoxs split with no validation/test samples)
    if len(dataloader) == 0:
        if flag:
            return 0.0, -999.0, -999.0, -999.0, -999.0, None
        else:
            return 0.0, -999.0, -999.0, -999.0, -999.0, None
    total_loss = 0.0
    total_loss1 = 0.0
    total_loss2 = 0.0
    total_loss3 = 0.0
    total_loss4 = 0.0
    num_batches = 0
    trt_raw_data_ls = []
    output2_ls = []
    ctl_raw_data_ls = []
    deg_pred_ls = []
    
    a, b, c, d = args.loss_weight
    MultiClassLoss = torch.nn.CrossEntropyLoss(reduction='sum')

    with torch.no_grad():
        for data in tqdm(dataloader, desc=f'Validate', leave=False):
            if scaler:
                with autocast(device_type='cuda'):
                    trt_output, ctl_output, deg_output, trt_raw_data, ctl_raw_data, _,cell_class_true, cell_class_predict = model(data)
                    num_samples =  trt_raw_data.shape[0]
                    num_batches += num_samples
                    loss1 = mse_loss_ls_sum(trt_output, trt_raw_data)
                    if cell_class_predict is not None:
                        cell_class_predict_1, cell_class_predict_2 = cell_class_predict
                        loss2 = MultiClassLoss(cell_class_predict_1, cell_class_true) + MultiClassLoss(cell_class_predict_2, cell_class_true)
                    else:
                        loss2 = mse_loss_ls_sum(ctl_output, ctl_raw_data)
                    loss3 = mse_loss_ls_sum(deg_output, trt_raw_data-ctl_raw_data)
                    loss4 = pcc_loss_sum(deg_output, trt_raw_data-ctl_raw_data)
                    weighted_loss = loss1 * a + loss2 * b + loss3 * c + loss4 * d
                    total_loss += weighted_loss.item()
                    total_loss1 += loss1.item()
                    total_loss2 += loss2.item()
                    total_loss3 += loss3.item()
                    total_loss4 += loss4.item()
            else:
                trt_output, ctl_output, deg_output, trt_raw_data, ctl_raw_data, _, cell_class_true, cell_class_predict = model(data)
                num_samples =  trt_raw_data.shape[0]
                num_batches += num_samples
                loss1 = mse_loss_ls_sum(trt_output, trt_raw_data)
                loss2 = mse_loss_ls_sum(ctl_output, ctl_raw_data)
                loss3 = mse_loss_ls_sum(deg_output, trt_raw_data-ctl_raw_data)
                loss4 = pcc_loss_sum(deg_output, trt_raw_data-ctl_raw_data)
                weighted_loss = loss1 * a + loss2 * b + loss3 * c + loss4 * d
                total_loss += weighted_loss.item()
                total_loss1 += loss1.item()
                total_loss2 += loss2.item()
                total_loss3 += loss3.item()
                total_loss4 += loss4.item()                        

            if flag:
                trt_raw_data_ls.append(trt_raw_data)
                output2_ls.append(deg_output + ctl_raw_data)
                ctl_raw_data_ls.append(ctl_raw_data)
                deg_pred_ls.append(deg_output)


    if flag:
        y_true = torch.cat(trt_raw_data_ls, dim=0).cpu().detach().numpy()
        y_pred = torch.cat(output2_ls, dim=0).cpu().detach().numpy()
        ctl_true = torch.cat(ctl_raw_data_ls, dim=0).cpu().detach().numpy()
        deg_pred = torch.cat(deg_pred_ls, dim=0).cpu().detach().numpy()
        metrics = get_metrics_new(y_true, y_pred, ctl_true)
        if args.output_profile:
            print('Saving output_profile....')
            output_profile_dict = {}
            output_profile_dict['y_true'] = y_true
            output_profile_dict['y_pred'] = y_pred
            output_profile_dict['ctl_true'] = ctl_true
            output_profile_dict['deg_true'] = y_true-ctl_true
            output_profile_dict['deg_pred'] = y_pred-ctl_true
            np.save(f'{expt_folder}/{fold}_predict_profile.npy', output_profile_dict)
    else:
        metrics = None

    if num_batches == 0:
        return 0.0, -999.0, -999.0, -999.0, -999.0, None

    avg_total_loss = total_loss / num_batches
    avg_loss1 = np.sqrt(total_loss1 / num_batches) if total_loss1 else -999
    if cell_class_predict is not None:
        avg_loss2 = (total_loss2 / num_batches) if total_loss2 else -999
    else:
        avg_loss2 = np.sqrt(total_loss2 / num_batches) if total_loss2 else -999
    avg_loss3 = np.sqrt(total_loss3 / num_batches) if total_loss3 else -999
    avg_loss4 = total_loss4 / num_batches if total_loss4 else -999



    return avg_total_loss, avg_loss1, avg_loss2, avg_loss3, avg_loss4, metrics





def infer(model, dataloader, args, config, flag=True, output_profile=False, scaler=None):
    model.eval()  
    trt_raw_data_ls = []
    output2_ls = []
    ctl_raw_data_ls = []
    pert_id_ls = []
    cell_id_ls = []
    pert_dose_ls = []
    pert_time_ls = []
    trt_attention_dict_ls= {}
    ctl_attention_dict_ls= {}
    cls_embed_ls = []


    with torch.no_grad():

        for data in tqdm(dataloader, desc=f'Infer', leave=False):
            trt_raw_data, ctl_raw_data, trt_raw_data_binned, ctl_raw_data_binned, drug_feat, pert_dose_idx, pert_time_idx, drug_idx, cell_idx, tissue_idx, pert_id, cell_id, pert_dose, pert_time = data

            if scaler:
                with autocast(device_type='cuda'):
                    trt_output, ctl_output, output_deg, _, _, attention_dict, _, _, cls_embed= model((trt_raw_data, ctl_raw_data, trt_raw_data_binned, ctl_raw_data_binned,drug_feat, pert_dose_idx, pert_time_idx, drug_idx, cell_idx, tissue_idx))
            else:
                trt_output, ctl_output, output_deg, _, _, attention_dict, _, _, cls_embed = model((trt_raw_data, ctl_raw_data, trt_raw_data_binned, ctl_raw_data_binned,drug_feat, pert_dose_idx, pert_time_idx, drug_idx, cell_idx, tissue_idx))                     

            if flag:
                trt_raw_data_ls.append(trt_raw_data)
                output2_ls.append(output_deg.cpu()+ctl_raw_data)
                ctl_raw_data_ls.append(ctl_raw_data)
                if args.output_attention:
                    trt_attention_dict, ctl_attention_dict = attention_dict
                    for key, value in trt_attention_dict.items():
                        if key not in trt_attention_dict_ls:
                            trt_attention_dict_ls[key] = []
                        if value is None:
                            print(f'value of {key} is None!')
                        trt_attention_dict_ls[key].append(value)
                    for key, value in ctl_attention_dict.items():
                        if key not in ctl_attention_dict_ls:
                            ctl_attention_dict_ls[key] = []
                        if value is None:
                            print(f'value of {key} is None!')
                        ctl_attention_dict_ls[key].append(value) 
                if args.output_cls_embed:
                    cls_embed_ls.append(cls_embed)


                pert_id = list(drug_idx)
                cell_id = list(cell_id)
                pert_dose = list(pert_dose)
                pert_time = list(pert_time)

                pert_id_ls.extend(pert_id)
                cell_id_ls.extend(cell_id)
                pert_dose_ls.extend(pert_dose)
                pert_time_ls.extend(pert_time)


        if flag:
            y_true = torch.cat(trt_raw_data_ls, dim=0).cpu().detach().numpy()
            y_pred = torch.cat(output2_ls, dim=0).cpu().detach().numpy()
            ctl_true = torch.cat(ctl_raw_data_ls, dim=0).cpu().detach().numpy()
            metrics, metrics_all_ls = get_metrics_infer(y_true, y_pred, ctl_true)

            metrics_all_ls['pert_id'] = np.array(pert_id_ls)
            metrics_all_ls['cell_id'] = np.array(cell_id_ls)
            metrics_all_ls['pert_dose'] = np.array(pert_dose_ls)
            metrics_all_ls['pert_time'] = np.array(pert_time_ls)
        else:
            metrics = None
            metrics_all_ls = None

        all_attention_dict = None
        if args.output_attention:
            all_attention_dict = {}
            all_trt_attention_dict = {}
            all_ctl_attention_dict = {}
            for key, value in trt_attention_dict_ls.items():
                all_trt_attention_dict[key] = torch.cat(value, dim=0).cpu().detach().numpy()
            for key, value in ctl_attention_dict_ls.items():
                all_ctl_attention_dict[key] = torch.cat(value, dim=0).cpu().detach().numpy()
            all_attention_dict['trt'] = all_trt_attention_dict
            all_attention_dict['ctl'] = all_ctl_attention_dict 

        all_cls_embed = None
        if args.output_cls_embed:
            all_cls_embed = torch.cat(cls_embed_ls, dim=0).cpu().detach().numpy()

        
        if output_profile:
            output_pred_profile = {}
            output_pred_profile['y_true'] = y_true
            output_pred_profile['y_pred'] = y_pred
            output_pred_profile['ctl_true'] = ctl_true
            output_pred_profile['deg_true'] = y_true-ctl_true
            output_pred_profile['deg_pred'] = y_pred-ctl_true
        else:
            output_pred_profile = None

    return metrics, metrics_all_ls, output_pred_profile, all_attention_dict, all_cls_embed






def main():

    args = arg_parse()

    # parse comma-separated args
    args.loss_weight = [float(x) for x in args.loss_weight.split(',')]

    start_time = time.time()
    timestamp = time.strftime('%Y%m%d_%H%M%S', time.localtime())

    config = load_config(args.config)

    if args.mode == 'infer':
        stage = 'infer'
    elif args.pretrained_model:
        stage = 'finetune'
    else:
        stage = 'pretrain'

    if stage == 'finetune' and not args.pretrain_dataset:
        try:
            ckpt = torch.load(args.pretrained_model, map_location='cpu')
            args.pretrain_dataset = ckpt.get('dataset', '')
        except Exception:
            pass

    if stage == 'finetune':
        dataset_id = f'{args.pretrain_dataset}-{args.dataset}'
    else:
        dataset_id = args.dataset

    expt_folder = str(get_config_root(args.config) / stage / 'output' / f'{args.dataset}_{args.split_method}')
    log_dir = os.path.join(expt_folder, 'logs')

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(os.path.join(log_dir, f'{timestamp}.log'))
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    logger.info('\n---------args-----------')
    log_nested_dict(vars(args), logger)
    logger.info('\n')

    set_random_seed(int(args.seed))

    logger.info('\n--------configs----------')
    log_nested_dict(config, logger)

    if torch.cuda.is_available():
        device = torch.device(args.device)
        torch.backends.cudnn.benchmark = True
        logger.info("Use GPU: %s\n" % device)
    else:
        device = torch.device("cpu")
        logger.info("Use CPU")
    

    if args.mode == 'train':
        logger.info('Train mode:')

        folds = config['split_methods'][args.split_method]['folds']

        for fold_idx, fold_name in enumerate(folds, start=1):

            tr_dataloader, val_dataloader, test_dataloader, adata = load_dataloader(args,config,logger,nfold=fold_name,return_rawdata=True)

            loss_func = None
            kl_loss = None

            if args.model =='XPert':
                model = XPertNet(args, config, device, logger)
            else:
                assert False, 'Model not found!'
            model.init_weights()
            model.to(device)
            logger.info(model)

            total_params = sum([param.nelement() for param in model.parameters()])
            logger.info("Number of parameters: %.2fM" % (total_params/1e6))

            if loss_func is not None:
                optimizer = torch.optim.Adam(
                                [
                                    {"params": loss_func.weights_param, "lr": 1e-7},
                                    {"params": model.parameters(), "lr": args.train_lr, "weight_decay": args.weight_decay}
                                ]
                            )
            else:
                optimizer =  torch.optim.Adam(model.parameters(), lr=args.train_lr, weight_decay=args.weight_decay)


            if args.lr_scheduler:
                lr_scheduler = torch.optim.lr_scheduler.LambdaLR(
                    optimizer, 
                    lr_lambda=lambda epoch: 1.0 if epoch < 40 else 0.5
                )
            else:
                lr_scheduler = None

            if args.use_gradscaler:
                gradscaler = GradScaler(device='cuda')
                logger.info('Using GradScaler to accelerate trainning..')
            else:
                gradscaler = None


            start_epoch = 0
            if args.resume_from:
                resume_path=args.resume_from
                logger.info(f'resume_from path:{resume_path}')
                checkpoint = torch.load(resume_path, map_location=device)
                model_dict = model.state_dict()
                pretrained_dict = {k:v for k,v in checkpoint['model_state_dict'].items() if k in model_dict}
                model_dict.update(pretrained_dict)
                model.load_state_dict(model_dict)
                start_epoch = checkpoint['epoch'] + 1
                logger.info(f'Load previous-trained parameters sucessfully! From epoch {start_epoch} to train……')
            if args.pretrained_model:
                pretrained_model_path=args.pretrained_model
                logger.info(f'pretrained model path:{pretrained_model_path}')
                checkpoint = torch.load(pretrained_model_path, map_location=device)
                model_dict = model.state_dict()

                exclude_keys = [k.strip() for k in args.exclude_keys.split(',') if k.strip()]
                pretrained_dict = {
                    k: v for k, v in checkpoint['model_state_dict'].items()
                    if k in model_dict
                    and v.shape == model_dict[k].shape
                    and not any(exclude_key in k for exclude_key in exclude_keys)
                }

                model_dict.update(pretrained_dict)
                model.load_state_dict(model_dict)
                start_epoch = 0
                logger.info(f'Load previous-trained parameters sucessfully! From epoch {start_epoch} to train……')

                finetune_params = args.finetune_params
                logger.info(f'Updating params: {finetune_params}')
                set_requires_grad(model, finetune_params=finetune_params)
                logger.info("Model Parameters:")
                for name, param in model.named_parameters():
                    logger.info(f"{name}: requires_grad={param.requires_grad}")

                optimizer = torch.optim.Adam(model.parameters(), lr=args.finetune_lr, weight_decay=args.finetune_weight_decay)

            
            max_epoch = args.num_epochs
            tr_loss_list = []
            val_loss_list = []
            model_name = f'{args.dataset}_{args.split_method}_fold-{fold_idx}'
            stopper = EarlyStopping(mode='lower', metric='mse', patience=args.patience, n_fold=fold_name, folder=expt_folder, dataset=model_name, logger=logger)

            early_stop = False
            for epoch in range(start_epoch, max_epoch):
                if epoch == start_epoch:
                    epoch_start_time = time.time()
                train_loss, train_loss1, train_loss2, train_loss3, train_loss4 = train(model, optimizer, tr_dataloader, args, config, scaler=gradscaler, epoch=epoch, loss_func=loss_func, kl_loss=kl_loss)
                tr_loss_list.append(train_loss)
                logger.info(f'Epoch {epoch}, Train Total Loss: {train_loss:.3f}, Loss1: {train_loss1:.3f}, Loss2: {train_loss2:.3f}, Loss3: {train_loss3:.3f}, Loss4: {train_loss4:.3f}')                       
                if epoch == start_epoch:    
                    logger.info(f'Epoch_{epoch}_train_time:{int(time.time() - epoch_start_time)/60:.2f}min')
                    valid_start_time = time.time()
                val_loss, val_loss1, val_loss2, val_loss3, val_loss4, _ = validate(model, val_dataloader, args, config, flag=False, scaler=gradscaler)
                early_stop = stopper.step(val_loss4, model, epoch, optimizer)
                val_loss_list.append(val_loss)
                logger.info(f'Epoch {epoch}, Valid Total Loss: {val_loss:.3f}, Loss1: {val_loss1:.3f}, Loss2: {val_loss2:.3f}, Loss3: {val_loss3:.3f}, Loss4: {val_loss4:.3f}')
                if early_stop or epoch == max_epoch-1:
                    best_model = stopper.load_checkpoint(model, optimizer)
                    best_epoch = epoch
                    break
                if lr_scheduler:
                    lr_scheduler.step()
                if epoch == start_epoch or epoch % 10 == 0:
                    current_lr = optimizer.param_groups[-1]['lr']
                    logger.info(f'Epoch {epoch+1}, Learning Rate: {current_lr}') 
                if epoch == start_epoch:
                    logger.info(f'Epoch_{epoch}_valid_time:{int(time.time() - valid_start_time)/60:.2f}min')


            try:
                save_loss_fig(tr_loss_list,val_loss_list,model_name,expt_folder,logger)
            except Exception as e:
                logger.info(f"An error occurred while saving the loss fig: {e}")

            best_model_path = os.path.join(expt_folder, f'{model_name}_best.pth')
            torch.save({
                'epoch': best_epoch,
                'model_state_dict': best_model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'best_val_loss': stopper.best_score,
                'dataset': args.dataset,
            }, best_model_path)
            logger.info(f'Best model saved to {best_model_path}')

            # Save best-val-loss metadata alongside model
            best_json_path = os.path.join(expt_folder, f'{model_name}_best.json')
            with open(best_json_path, 'w') as f:
                json.dump({'best_val_loss': stopper.best_score, 'timestamp': timestamp}, f)

            logger.info(f'{args.split_method}_fold-{fold_idx} training is done! Training_time:{int(time.time() - start_time)/60:.2f}min')
            logger.info('Start testing ... ')

            logger.info('Tesing on train dataset:')
            _,_,_,_,_, tr_metrics = validate(best_model, tr_dataloader, args, config, flag=True, scaler=gradscaler)
            logger.info('Training Metrics:')
            if tr_metrics is not None:
                log_nested_dict(tr_metrics, logger, indent=1)

            logger.info('Tesing on Valid dataset:')
            _,_,_,_,_, val_metrics = validate(best_model, val_dataloader, args, config, flag=True, scaler=gradscaler)
            logger.info('Validation Metrics:')
            if val_metrics is not None:
                log_nested_dict(val_metrics, logger, indent=1)

            logger.info('Tesing on test dataset:')
            _,_,_,_,_, test_metrics = validate(model, test_dataloader, args, config, flag=True, scaler=gradscaler)
            test_metrics = val_metrics
            logger.info('Test Metrics:')
            if test_metrics is not None:
                log_nested_dict(test_metrics, logger, indent=1)

            logger.info('Saving result……')
            result_dict_list = [m for m in [tr_metrics, val_metrics, test_metrics] if m is not None]
            result_df = pd.DataFrame(result_dict_list) if result_dict_list else pd.DataFrame()

            try:
                result_df.to_csv(os.path.join(expt_folder, f'{model_name}_train_result.csv'), index=False)
                logger.info('Results successfully saved!!!\n\n')
            except Exception as e:
                logger.info(f"An error occurred while saving the [{model_name}_train_result.csv] file: {e}")

        logger.info(f'All folds training are done! Training time: {int(time.time() - start_time)/60:.2f} min')
            


    elif args.mode == 'test':
        logger.info('Test mode...')
        
        folds = config['split_methods'][args.split_method]['folds']

        for fold_idx, fold_name in enumerate(folds, start=1):
            test_dataloader = load_test_dataloader(args, config, logger, fold_name)

            model_name = f'{args.dataset}_{args.split_method}_fold-{fold_idx}'

            if args.model =='XPert':
                model = XPertNet(args, config, device, logger)
            else:
                assert False, 'Model not found!'
            model.init_weights()
            model.to(device)

            if args.saved_model:
                saved_model = args.saved_model
            elif args.saved_model_path:
                saved_model = osp.join(args.saved_model_path, f'{model_name}_best.pth')
            else:
                logger.info('No trained parameters are provided!')
            checkpoint = torch.load(saved_model)
            model.load_state_dict(checkpoint['model_state_dict'])

            if args.use_gradscaler:
                gradscaler = GradScaler(device='cuda')
                logger.info('Using GradScaler to accelerate trainning..')
            else:
                gradscaler = None


            logger.info('Start testing ... ')
            start_test_time = time.time()
            logger.info('Tesing on test dataset:')
            _,_,_,_,_, test_metrics = validate(model, test_dataloader, args, config, flag=True, scaler=gradscaler, expt_folder=expt_folder, fold=fold_name)
            logger.info('Test Metrics:')
            log_nested_dict(test_metrics, logger, indent=1)

            logger.info(f'Test time: {int(time.time() - start_test_time)/60:.2f} min')

            test_metrics = test_metrics
            val_metrics = test_metrics
            tr_metrics = test_metrics

            logger.info('Saving result……')
            result_dict_list = [tr_metrics, val_metrics, test_metrics]
            result_df = pd.DataFrame(result_dict_list)

            try:
                result_df.to_csv(os.path.join(expt_folder, f'{model_name}_train_result.csv'), index=False)
                logger.info('Results successfully saved!!!\n\n')
            except Exception as e:
                logger.info(f"An error occurred while saving the [{model_name}_train_result.csv] file: {e}")

    elif args.mode == 'infer':

        logger.info('Infer mode...')

        folds = config['split_methods'][args.split_method]['folds']

        for fold_idx, fold_name in enumerate(folds, start=1):
            test_dataloader = load_infer_dataloader(args, config, logger, fold_name)

            model_name = f'{args.dataset}_{args.split_method}_fold-{fold_idx}'

            if args.model =='XPert':
                model = XPertNet(args, config, device, logger)
            else:
                assert False, 'Model not found!'
            
            model.to(device)

            total_params = sum([param.nelement() for param in model.parameters()])
            logger.info("Number of parameters: %.2fM" % (total_params/1e6))

            if args.saved_model:
                saved_model = args.saved_model
            elif args.saved_model_path:
                saved_model = osp.join(args.saved_model_path, f'{model_name}_best.pth')
            else:
                logger.info('No trained parameters are provided!')

            checkpoint = torch.load(saved_model)

            model_dict = model.state_dict()
            pretrained_dict = {k:v for k,v in checkpoint['model_state_dict'].items() if k in model_dict and v.shape == model_dict[k].shape}
            model_dict.update(pretrained_dict)
            model.load_state_dict(model_dict)


            logger.info('Start testing ... ')
           
            logger.info('Tesing on test dataset:')
            test_metrics, test_metrics_all_ls, test_output_pred_profile, test_all_attention_dict, all_cls_embed = infer(model, test_dataloader, args, config, flag=True, output_profile=args.output_profile, scaler=args.use_gradscaler)

            logger.info('Test Metrics:')
            log_nested_dict(test_metrics, logger, indent=1) 

                     
            logger.info('Saving result……')

            result_df = pd.DataFrame.from_dict(test_metrics, orient='index').T         
            try:
                result_df.to_csv(expt_folder + f'/{model_name}_train_result.csv')
                logger.info('Results successfully saved!!!\n\n')
            except Exception as e:
                logger.info(f"An error occurred while saving the [{model_name}_train_result.csv] file: {e}")

            if test_metrics_all_ls is not None:
                result_ls_df = pd.DataFrame.from_dict(test_metrics_all_ls)          
                try:
                    result_ls_df.to_csv( expt_folder + f'/{model_name}_test_samples_result.csv' )
                    logger.info('Results successfully saved!!!\n\n')
                except Exception as e:
                    logger.info(f"An error occurred while saving the [{model_name}_test_result_all_samples.csv] file: {e}")

            if test_output_pred_profile is not None:
                try:
                    np.save( expt_folder + f'/{model_name}_test_samples_prediction_profile.npy' , test_output_pred_profile)
                    logger.info('Profiles successfully saved!!!\n\n')
                except Exception as e:
                    logger.info(f"An error occurred while saving the [{model_name}_test_samples_prediction_profile.npy] file: {e}")

            if args.output_attention:
                try:
                    output_file_path = expt_folder + f'/{model_name}_test_samples_attention_scores.pkl'
                    with open(output_file_path, "wb") as f:
                        pickle.dump(test_all_attention_dict, f)
                    logger.info('Attention scores successfully saved!!!')
                except Exception as e:
                    logger.info(f"An error occurred while saving the [{model_name}_test_samples_attention_scores.npy] file: {e}")

            if args.output_cls_embed:
                try:
                    np.save( expt_folder + f'/{model_name}_test_samples_cls_embed.npy' , all_cls_embed)
                    logger.info('Profiles successfully saved!!!\n\n')
                except Exception as e:
                    logger.info(f"An error occurred while saving the [{model_name}_test_samples_cls_embed.npy] file: {e}") 


if __name__ == '__main__':
   
    main()
