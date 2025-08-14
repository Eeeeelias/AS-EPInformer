# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import pandas as pd
# torch
import torch
from torch import nn
import torch.utils.data as data_utils
from torch.utils.data import Subset
from torch.nn.utils import parameters_to_vector
from torch.nn.functional import cosine_similarity
from torchviz import make_dot

from denseweight import DenseWeight

from scipy import stats
from sklearn.metrics import matthews_corrcoef, mean_squared_error, r2_score, roc_auc_score, balanced_accuracy_score
from sklearn.model_selection import train_test_split
from tqdm import tqdm
#
#from softadapt import LossWeightedSoftAdapt


# logging
# from model.EPInformer import EPInformer_v2, enhancer_predictor_256bp
from scripts.EarlyStopping import EarlyStopping, EarlyStoppingMulti
from scripts.loss_functions import dense_loss

def get_lr(optimizer):
    for param_group in optimizer.param_groups:
        return param_group['lr']


def combine_hurdle_outputs(binary_logits, splicing_pred, threshold=0.5):
    """Combine binary and regression outputs from a hurdle model."""
    binary_probs = torch.sigmoid(binary_logits)
    is_not_1 = (binary_probs > threshold).squeeze(-1)  # True = not 1, False = 1

    final_pred = torch.where(is_not_1, splicing_pred.squeeze(-1), torch.ones_like(splicing_pred.squeeze(-1)))

    return list(final_pred.cpu().detach().numpy())


def get_sample_weights(trainloader, device='cpu', filter_ones=True):
    dw = DenseWeight(alpha=0.6)
    psi_reg = []
    psi_binary = []
    for data in trainloader:
        _, _, _, _, _, _, y_psi, _ = data
        flat_psi = y_psi.flatten().cpu().numpy()
        # add to psi binary 1 if the value is 1.0, else 0
        psi_binary.extend((flat_psi == 1.0).astype(int))
        # filter out 1.0 values since they are protected by the hurdle
        if filter_ones:
            flat_psi = flat_psi[flat_psi != 1.0]
            if len(flat_psi) > 0:
                psi_reg.extend(flat_psi)
        else:
            psi_reg.extend(flat_psi)
    dw.fit(np.array(psi_reg))

    num_positive = np.sum(np.array(psi_binary) == 1)
    num_negative = len(psi_binary) - num_positive
    pos_weight = torch.tensor(num_negative / num_positive, dtype=torch.float32, device=device)
    return dw, pos_weight


def get_sample_weights_binary(trainloader, device='cpu'):
    psi_binary = []
    for data in trainloader:
        _, y_psi, _ = data
        flat_psi = y_psi.flatten().cpu().numpy()
        mask = (flat_psi <= 0.2) | (flat_psi >= 0.8)
        psi_masked = flat_psi[mask]
        psi_binary.extend((psi_masked > 0.5).astype(int))
    num_positive = np.sum(np.array(psi_binary) == 1)
    num_negative = len(psi_binary) - num_positive
    pos_weight = torch.tensor(num_negative / num_positive, dtype=torch.float32, device=device)
    print(f"Number of positive samples: {num_positive}, Number of negative samples: {num_negative}\n" \
        f"Positive weight: {pos_weight}")

    return pos_weight


def get_loss_function(expr_loss_type='mse', splice_loss_type='bce', **kwargs):
    """Get the loss function based on the specified type."""
    if expr_loss_type == 'mse':
        L_expr = nn.MSELoss()
    elif expr_loss_type == 'smoothl1':
        L_expr = nn.SmoothL1Loss()
    else:
        raise ValueError(f"Unsupported expression loss type: {expr_loss_type}")

    if splice_loss_type == 'bce':
        L_splice = nn.BCEWithLogitsLoss(**kwargs)
    elif splice_loss_type == 'smoothl1':
        L_splice = nn.SmoothL1Loss(**kwargs)
    else:
        raise ValueError(f"Unsupported splice loss type: {splice_loss_type}")
    return L_expr, L_splice


def combine_losses(loss_expr, loss_splice, predict_type='multi', weights=(1.0, 1.0)) -> torch.Tensor:
    loss = torch.tensor(0.0, device=loss_expr.device)
    if predict_type == 'multi':
        if weights[0] == -np.inf: # expression loss should not be used
            loss = loss_splice
        if weights[1] == -np.inf: # splice loss should not be used
            loss = loss_expr
        if not (weights[0] == -np.inf or weights[1] == -np.inf): # both losses should be used
            # if both weights are not -inf, combine the losses
            loss = weights[0] * loss_expr + weights[1] * loss_splice

    elif predict_type == 'splice':
        loss = loss_splice

    elif predict_type == 'RNA':
        loss = loss_expr
    else:
        raise ValueError(f"Unsupported predict type: {predict_type}")
        
    return loss


def setup_loss_file(saved_model_path):
    if saved_model_path and not os.path.exists(saved_model_path):
        os.mkdir(saved_model_path)
    if saved_model_path and not os.path.exists(saved_model_path + "/losses.csv"):
        loss_file = open(saved_model_path + "/losses.csv", "w", encoding='utf-8')
        loss_file.write("fold,epoch,training_loss,train_expr_loss,train_splice_loss,val_mse,val_r2,val_loss,val_expr,val_splice\n")
    elif saved_model_path:
        loss_file = open(saved_model_path + "/losses.csv", "a", encoding='utf-8')
    else:
        loss_file = None
    return loss_file


def binarize_psi(pred_psi, true_psi):
    # filter out samples where true_psi is > 0.2 and < 0.8
    mask = (true_psi <= 0.2) | (true_psi >= 0.8)

    bin_pred = pred_psi[mask]
    bin_true = true_psi[mask]

    bin_true = (bin_true >= 0.8).long().float()

    # return 0, 0 if it happens that there are no samples left after filtering
    if len(bin_pred) == 0 or len(bin_true) == 0:
        print("No samples left after filtering, returning tensors of zeros")
        return torch.tensor([]), torch.tensor([])

    return bin_pred, bin_true



def train(net, training_dataset, fold_i, saved_model_path='../models', learning_rate=1e-4, model_logger=None, 
          fixed_encoder = False, n_enhancers = 50, valid_dataset = None, model_name = '', batch_size = 64, 
          device = 'cuda', stratify=None, epochs=100, valid_size=1000, predict='multi', loss_class=None, 
          weigh_samples=False, expr_loss_type='mse', splice_loss_type='bce'):
    
    loss_file = setup_loss_file(saved_model_path)

    train_ds = training_dataset
    valid_ds = valid_dataset

    # fix encoder parameter
    if fixed_encoder:
        print('fixed parameter of encoder')
        for name, value in net.named_parameters():
            if name.startswith('seq_encoder'):
                value.requires_grad = False

    print("fold", fold_i ,"training data:", len(train_ds), "validated data:", len(valid_ds), 'total data:', len(training_dataset))
    trainloader = data_utils.DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=5, pin_memory=True)
    if saved_model_path is not None:
        early_stopping = EarlyStoppingMulti(expr_patience=1, splice_patience=20,
                                            path=f"{saved_model_path}/fold_{fold_i}_best_{model_name}_checkpoint.pt",
                                            verbose=True)
    else:
        early_stopping = EarlyStoppingMulti(expr_patience=3, splice_patience=3, verbose=True)

    
    # get all PSI values from training dataset
    if predict == 'multi' and weigh_samples:
        dw, pos_weight = get_sample_weights(trainloader, device=device)
        print('Dense weight:', dw)
        print('Positive weight:', pos_weight)
    else:
        dw = None
        pos_weight = None

    # Loss functions
    reduction = 'mean' if not weigh_samples else 'none'
    # L_expr, L_splice = get_loss_function(expr_loss_type=expr_loss_type, splice_loss_type=splice_loss_type, reduction=reduction)
    # BINARY TEST ONLY
    pos_weight = get_sample_weights_binary(trainloader=trainloader, device=device)
    L_binary = nn.BCEWithLogitsLoss(reduction=reduction, pos_weight=pos_weight)
    learned_loss = True if loss_class is not None else False
    loss_weights = [0.5, 1.0] # modified by early stopping

    # optimizer
    all_params = net.parameters() if not learned_loss else list(net.parameters()) + list(loss_class.parameters()) # type: ignore
    optimizer = torch.optim.AdamW(all_params, lr=learning_rate, weight_decay=1e-6)
    net.train()

    # when using the xpu device, optimise the model for the xpu
    if device == 'xpu':
        import intel_extension_for_pytorch as ipex
        print("Using XPU device for training")
        net, optimizer = ipex.optimize(net, optimizer=optimizer)

    print('Model name:', net.name)
    lrs = []

    # softadapt init stuff
    # soft_adapt = LossWeightedSoftAdapt(beta=0.1)
    # update_epochs = 2
    losses_expr = []
    losses_splice = []
    # adapt_weights = torch.tensor([1.0, 1.0], device=device)

    for epoch in range(epochs):
        net.train()
        print('learning rate:', get_lr(optimizer))
        running_loss = 0
        expression_loss = 0
        splice_loss = 0
        #if epoch % update_epochs == 0 and epoch > 0:
        #    # save both tensors
        #    adapt_weights = soft_adapt.get_component_weights(torch.tensor(losses_expr), torch.tensor(losses_splice))  
        # print(f"Expr Adapt Weight: {adapt_weights[0]:.4f}, Splice Adapt Weight: {adapt_weights[1]:.4f}")

        # print('model training mode is:', net.training)
        for data in tqdm(trainloader, bar_format='{desc:<8}{percentage:3.0f}%|{bar:20}{r_bar}', desc=f'Epoch {epoch + 1}'):
            # print(inputs.size())
            optimizer.zero_grad()
            input_ex, y_psi, _ = data
            input_ex = input_ex.float().to(device)

            y_psi = y_psi.float().to(device)

            # pred_expr, pred_splice_binary, pred_splice, _ = net(input_pe, input_ex, p_input_hist, input_pe_feat, input_cell)
            pred_splice_binary = net(input_ex)

            # loss_expr = L_expr(pred_expr, y_expr.reshape(pred_expr.shape))
            # loss_splice = L_splice(pred_splice, y_psi.reshape(pred_splice.shape))

            pred_bin, y_bin = binarize_psi(pred_splice_binary, y_psi)
            if len(pred_bin) == 0 or len(y_bin) == 0:
                continue
            loss_binary = L_binary(pred_bin, y_bin.reshape(pred_bin.shape))

            #if weigh_samples:
            #    loss_splice = dense_loss(pred_splice, y_psi.reshape(pred_splice.shape), dw, loss_fn=L_splice)

            # loss = combine_losses(loss_expr, loss_splice, predict_type=predict, weights=loss_weights)
            loss = loss_binary # BINARY TEST ONLY
            if not os.path.exists("images/binary_graph.png"):
                print("Creating computation graph for the first time")
                make_dot(loss, params=dict(net.named_parameters())).render("images/binary_graph", format="png")
            # propagate the loss backward
            loss.backward()
            # update the gradients
            optimizer.step()
            running_loss += loss.item()

            # expression_loss += loss_expr.item()
            # splice_loss += loss_splice.item()
            expression_loss += 0 # CHANGED BECAUSE OF BINARY
            splice_loss += 0 # CHANGED BECAUSE OF BINARY

        print(f"[Epoch {epoch + 1}] overall loss: {running_loss/len(trainloader):.5f}, "
              f"expression loss: {expression_loss/len(trainloader):.5f}, " \
              f"splice loss: {splice_loss/len(trainloader):.5f}")

        # softadapt updates
        losses_expr.append(expression_loss/len(trainloader))
        losses_splice.append(splice_loss/len(trainloader))

        val = validate(net, valid_ds, n_enhancers=n_enhancers, device=device, predict=predict, 
                        expr_loss_type=expr_loss_type, splice_loss_type=splice_loss_type)

        print('Validation R square all:', val['r2'])
        early_stopping(val['expression_loss'], val['splice_loss'], net, epoch)
        if loss_file is not None:
            loss_file.write(f"{fold_i},{epoch+1},{running_loss/len(trainloader)},{expression_loss/len(trainloader)}," \
                            f"{splice_loss/len(trainloader)},{val['mse']},{val['r2']},{val['total_loss']}," \
                            f"{val['expression_loss']},{val['splice_loss']}\n")
            
            loss_file.flush()
        if model_logger is not None:
            label_type = net.name.split('.')[-1]
            model_logger.add([fold_i, epoch, running_loss/len(trainloader), val['mse'], val['pearsonr'], val['r2'], 
                              val['pearsonr'], val['r2'], early_stopping.counter, label_type])
            # model_logger.save("./EPInfomrer_log/{}.crossValid.log".format(net.name.replace('.'+label_type, '')))
        if early_stopping.expr_early_stop and early_stopping.splice_early_stop:
            print("Early stopping")
            break
        elif early_stopping.expr_early_stop:
            print("Early stopping on expression loss")
            loss_weights[0] = -np.inf
        elif early_stopping.splice_early_stop:
            print("Early stopping on splice loss")
            loss_weights[1] = -np.inf

    return lrs

def validate(net, valid_ds,  net_type = 'seq_feat_dist', n_enhancers=50, batch_size=16, device = 'cuda', 
             predict='multi', loss_weights=(1.0, 1.0, 1.0), expr_loss_type='mse', splice_loss_type='bce') -> dict:
    validloader = data_utils.DataLoader(valid_ds, batch_size=batch_size, pin_memory=True, num_workers=0)
    net.eval()
    L_expr, L_psi = get_loss_function(expr_loss_type=expr_loss_type, splice_loss_type=splice_loss_type)
    L_binary = nn.BCEWithLogitsLoss(reduction='mean')

    with torch.no_grad():
        preds = []
        actual = []
        preds_psi = []
        actual_psi = []
        loss_e = 0
        expression_loss = 0
        splice_loss = 0
        for data in validloader:
            # print(inputs.size())
            input_seg, y_psi, _ = data
            input_seg = input_seg.float().to(device)
            y_psi = y_psi.float().to(device)
            # print(input_P.shape, input_E.shape, input_Emask.shape)

            # pred_expr, pred_splice_binary, pred_splice, _ = net(input_pe, input_seg, input_feat, input_dist, input_cell)
            pred_splice_binary = net(input_seg)


            # outputs = list(pred_expr.flatten().cpu().detach().numpy())
            outputs = list(pred_splice_binary.flatten().cpu().detach().numpy()) # TEMP FOR SIMPLE BINARY
            labels = list(pred_splice_binary.flatten().cpu().detach().numpy())

            #outputs_psi = list(torch.sigmoid(pred_splice).flatten().cpu().detach().numpy())
            #labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            #loss_expr = L_expr(pred_expr, y_expr.reshape(pred_expr.shape))
            expression_loss += 0 # CHANGED BECAUSE OF BINARY
            
            #loss_splice = L_psi(pred_splice, y_psi.reshape(pred_splice.shape))
            # splice_loss += loss_splice.item()

            # BINARY TEST ONLY
            pred_bin, y_bin = binarize_psi(pred_splice_binary, y_psi)
            if len(pred_bin) == 0 or len(y_bin) == 0:
                continue
            loss_binary = L_binary(pred_bin, y_bin.reshape(pred_bin.shape))
            splice_loss += loss_binary.item()
            loss_e += loss_binary
            outputs_psi = list(torch.sigmoid(pred_bin).flatten().cpu().detach().numpy())
            labels_psi = list(y_bin.flatten().cpu().detach().numpy())

            # loss_e += combine_losses(loss_expr, loss_splice, predict_type=predict)

            preds += outputs
            actual += labels
            preds_psi += outputs_psi
            actual_psi += labels_psi

    min_max_expr = [np.min(preds), np.max(preds)]
    min_max_psi = [np.min(preds_psi), np.max(preds_psi)]

    try:
        r2_value = r2_score(actual, preds)
        peasonr, _ = stats.pearsonr(preds, actual)
        mse = mean_squared_error(actual, preds)
    except ValueError:
        r2_value = 0
        peasonr = 0
        mse = 0
    print(f"[Validat] overall loss: {loss_e / len(validloader):.5f}, "
          f"expression loss: {expression_loss / len(validloader):.5f}, " \
          f"splice loss: {splice_loss / len(validloader):.5f}")
    print("\n### Validation ### TPM expresion ###")
    print(f'### Min-Max Expression: {min_max_expr[0]:.5f}, {min_max_expr[1]:.5f}')
    print(f"### MSE:{mse:.5f} R²: {r2_value:.5f} PeasonR: {peasonr:.5f}\n")

    try:
        r2_value_psi = r2_score(actual_psi, preds_psi)
        peasonr_psi, _ = stats.pearsonr(preds_psi, actual_psi)
        mse_psi = mean_squared_error(actual_psi, preds_psi)
        # BINARY TEST ONLY
        r2_value_psi = 0
        peasonr_psi = 0
        mse_psi = 0
        auroc = roc_auc_score(actual_psi, preds_psi)
        mcc = matthews_corrcoef(actual_psi, (np.array(preds_psi) >= 0.5).astype(int))
        bal_acc = balanced_accuracy_score(actual_psi, (np.array(preds_psi) >= 0.5).astype(int))
    except ValueError as e:
        print(f"ValueError in PSI validation, setting metrics to 0: {e}")
        r2_value_psi, peasonr_psi, mse_psi = 0, 0, 0
        auroc, mcc, bal_acc = 0, 0, 0
    print("### Validation ### PSI expression ###")
    #print(f'### Min-Max PSI: {min_max_psi[0]:.5f}, {min_max_psi[1]:.5f}')
    #print(f"### MSE:, {mse_psi:.5f} R²: {r2_value_psi:.5f} PeasonR: {peasonr_psi:.5f}")
    print(f"### AUROC: {auroc:.5f} MCC: {mcc:.5f} Bal Acc: {bal_acc:.5f}\n")

    # get average r2 
    if predict == 'multi':
        avg_r2 = (r2_value + r2_value_psi) / 2
    elif predict == 'splice':
        avg_r2 = r2_value_psi
    else:
        avg_r2 = r2_value

    return {'mse': mse, 'r2': avg_r2, 'peasonr': peasonr, 'total_loss': loss_e / len(validloader), 
            'expression_loss': expression_loss / len(validloader), 'splice_loss': splice_loss / len(validloader)}


def test(net, test_ds, fold_i, model_name = None, saved_model_path=None, batch_size=64, device = 'cuda', 
         model_type='best', normals=None, predict='multi'):
    testloader = data_utils.DataLoader(test_ds, batch_size=batch_size, pin_memory=True, num_workers=0)

    if saved_model_path is not None:
        checkpoint = torch.load(saved_model_path + "/fold_" + str(fold_i) + "_best_"+model_name+"_checkpoint.pt", 
                                weights_only=False)
        net.load_state_dict(checkpoint['model_state_dict'])
        print(model_name,'loaded!')
        
    net.eval()
    with torch.no_grad():
        preds = []
        actual = []
        preds_psi = []
        actual_psi = []
        cell_line = []
        uncorr_pred_expr = []
        uncorr_pred_splice = []
        uncorr_actual_expr = []
        uncorr_actual_psi = []

        ensid_list = []
        for data in tqdm(testloader):
            input_seg, y_psi, eid = data
            input_seg = input_seg.float().to(device)
            y_psi = y_psi.float().to(device)
            #pred_expr, pred_splice_binary, pred_splice, _ = net(input_pe, input_seg, input_feat, input_dist, input_cell)
            pred_splice_binary = net(input_seg)
            pred_expr, pred_splice = pred_splice_binary, pred_splice_binary


            if normals:
                uncorr_pred_ep = pred_expr.clone()
                uncorr_pred_sp = pred_splice.clone()
                uncorr_y_expr = y_expr.clone()
                uncorr_y_psi = y_psi.clone()

                pred_expr = (pred_expr * normals['std_expr']) + normals['mean_expr']
                y_expr = (y_expr * normals['std_expr']) + normals['mean_expr']
          
                corr_pred_splice = []
                corr_y_psi = []
                for p_splice, y_splice in zip(pred_splice, y_psi):
                    if y_splice != 1.0:
                        p_splice = (p_splice * normals['std_psi']) + normals['mean_psi']
                        y_splice = (y_splice * normals['std_psi']) + normals['mean_psi']
                    corr_pred_splice.append(p_splice)
                    corr_y_psi.append(y_splice)

                pred_splice = torch.tensor(corr_pred_splice, device=device)
                y_psi = torch.tensor(corr_y_psi, device=device)

            outputs = list(pred_expr.flatten().cpu().detach().numpy())
            labels = list(pred_expr.flatten().cpu().detach().numpy())
            cell_lines = ['K562'] * len(outputs)

            outputs_psi = list(torch.sigmoid(pred_splice).flatten().cpu().detach().numpy())
            labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            # BINARY TEST ONLY
            outputs_psi = list(torch.sigmoid(pred_splice_binary).flatten().cpu().detach().numpy())
            labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            preds_psi += outputs_psi
            actual_psi += labels_psi

            preds += outputs
            actual += labels
            ensid_list += eid
            cell_line += cell_lines

            if normals:
                uncorr_pred_expr += list(uncorr_pred_ep.flatten().cpu().detach().numpy())
                uncorr_pred_splice += list(uncorr_pred_sp.flatten().cpu().detach().numpy())
                uncorr_actual_expr += list(uncorr_y_expr.flatten().cpu().detach().numpy())
                uncorr_actual_psi += list(uncorr_y_psi.flatten().cpu().detach().numpy())


    sys.stdout.flush()
    df = pd.DataFrame(index=np.array(ensid_list).flatten())
    df['PredExpr'] = preds
    df['ActualExpr'] = actual
    df['CellLine'] = cell_line
    if len(preds) == len(preds_psi):
        df['PredPsi'] = preds_psi
        df['ActualPsi'] = actual_psi
    if len(preds) == len(uncorr_pred_expr):
        df['UncorrPredExpr'] = uncorr_pred_expr
        df['UncorrActualExpr'] = uncorr_actual_expr
        df['UncorrPredPsi'] = uncorr_pred_splice
        df['UncorrActualPsi'] = uncorr_actual_psi
    df['fold_idx'] = fold_i

    if saved_model_path is not None:
        df.to_csv(saved_model_path + "/fold_" + str(fold_i) + "_"+ model_name + "_predictions.csv")
    return df
