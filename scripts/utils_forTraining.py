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
from denseweight import DenseWeight

from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from tqdm import tqdm
# logging
# from model.EPInformer import EPInformer_v2, enhancer_predictor_256bp
from scripts.loss_functions import hurdle_loss

def get_lr(optimizer):
    for param_group in optimizer.param_groups:
        return param_group['lr']


def combine_hurdle_outputs(binary_logits, splicing_pred, threshold=0.5):
    """Combine binary and regression outputs from a hurdle model."""
    binary_probs = torch.sigmoid(binary_logits)
    is_not_1 = (binary_probs > threshold).squeeze(-1)  # True = not 1, False = 1

    final_pred = torch.where(is_not_1, splicing_pred.squeeze(-1), torch.ones_like(splicing_pred.squeeze(-1)))

    return list(final_pred.cpu().detach().numpy())


def get_sample_weights(trainloader, device='cpu'):
    dw = DenseWeight(alpha=0.6)
    psi_reg = []
    psi_binary = []
    for data in trainloader:
        _, _, _, _, _, y_psi, _ = data
        flat_psi = y_psi.flatten().cpu().numpy()
        # add to psi binary 1 if the value is 1.0, else 0
        psi_binary.extend((flat_psi == 1.0).astype(int))
        # filter out 1.0 values since they are protected by the hurdle
        flat_psi = flat_psi[flat_psi != 1.0]
        if len(flat_psi) > 0:
            psi_reg.extend(flat_psi)
    dw.fit(np.array(psi_reg))

    num_positive = np.sum(np.array(psi_binary) == 1)
    num_negative = len(psi_binary) - num_positive
    pos_weight = torch.tensor(num_negative / num_positive, dtype=torch.float32, device=device)
    return dw, pos_weight


class Logger():
    """A logging class that can report or save metrics.

    This class contains a simple utility for saving statistics as they are
    generated, saving a report to a text file at the end, and optionally
    print the report to screen one line at a time as it is being generated.
    Must begin using the `start` method, which will reset the logger.

    Parameters
    ----------
    names: list or tuple
        An iterable containing the names of the columns to be logged.

    verbose: bool, optional
        Whether to print to screen during the logging.
    """

    def __init__(self, names, verbose=False):
        self.names = names
        self.verbose = verbose
        self.data = {name: [] for name in self.names}

    def start(self):
        """Begin the recording process."""

        if self.verbose:
            print("\t".join(self.names))

    def add(self, row):
        """Add a row to the log.

        This method will add one row to the log and, if verbosity is set,
        will print out the row to the log. The row must be the same length
        as the names given at instantiation.

        Parameters
        ----------
        args: tuple or list
            An iterable containing the statistics to be saved.
        """

        assert len(row) == len(self.names)

        for name, value in zip(self.names, row):
            self.data[name].append(value)

        if self.verbose:
            print("\t".join(map(str, [round(x, 4) if isinstance(x, float) else x
                for x in row])))

    def save(self, name):
        """Write a log to disk.


        Parameters
        ----------
        name: str
            The filename to save the logs to.
        """
        pd.DataFrame(self.data).to_csv(name, sep='\t', index=False)

class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, patience=3, verbose=False, delta=0, path='checkpoint.pt'):
        """
        Args:
            patience (int): How long to wait after last time validation loss improved.
                            Default: 6
            verbose (bool): If True, prints a message for each validation loss improvement.
                            Default: False
            delta (float): Minimum change in the monitored quantity to qualify as an improvement.
                            Default: 0
            path (str): Path for the checkpoint to be saved to.
                            Default: 'checkpoint.pt'
        """
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.inf
        self.delta = delta
        self.path = path

    def __call__(self, val_loss, model, epoch_i):
        score = -val_loss
        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model, epoch_i)
        elif score < self.best_score + self.delta:
            self.counter += 1
            print(f'EarlyStopping counter: {self.counter} out of {self.patience}', 'best_score', self.best_score)
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model, epoch_i)
            self.counter = 0

    def save_checkpoint(self, val_loss, model, epoch_i):
        '''Saves model when validation loss decrease.'''
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')

        torch.save({
                'epoch': epoch_i,
                'model_state_dict': model.state_dict(),
                'loss': val_loss,
                },
                self.path)
        print('Saving ckpt at', self.path)
        # torch.save(model.state_dict(), self.path)
        self.val_loss_min = val_loss


def train(net, training_dataset, fold_i, saved_model_path='../models', learning_rate=1e-4, model_logger=None, 
          fixed_encoder = False, n_enhancers = 50, valid_dataset = None, model_name = '', batch_size = 64, 
          device = 'cuda', stratify=None, epochs=100, valid_size=1000, predict='multi', loss_class=None, 
          weigh_samples=False):
    if not os.path.exists(saved_model_path):
        os.mkdir(saved_model_path)
    if not os.path.exists(saved_model_path + "/losses.csv"):
        loss_file = open(saved_model_path + "/losses.csv", "w", encoding='utf-8')
        loss_file.write("fold,epoch,training_loss,expresssion_loss,splice_loss,validation_mse,validation_r2,val_loss_total,val_loss_expr,val_loss_splice\n")
    else:
        loss_file = open(saved_model_path + "/losses.csv", "a", encoding='utf-8')

    if valid_dataset is not None:
        # genereate 1024 random indices for training dataset, i.e., from range(len(training_dataset))
        train_ds = training_dataset # Subset(training_dataset, random.sample(range(len(training_dataset)), 1024))
        valid_ds = valid_dataset
    else:
        train_idx, val_idx = train_test_split(list(range(len(training_dataset))), test_size=valid_size,
                                              shuffle=True, random_state=66, stratify=stratify)
        train_ds = Subset(training_dataset, train_idx)
        valid_ds = Subset(training_dataset, val_idx)
    # print quantiles of the training response and the validation response
    print('Training expression quantiles:', np.quantile([x[4].numpy() for x in train_ds], [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))
    print('Validation expression quantiles:', np.quantile([x[4].numpy() for x in valid_ds], [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))
    print('Training PSI quantiles:', np.quantile([x[5].numpy() for x in train_ds], [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))
    print('Validation PSI quantiles:', np.quantile([x[5].numpy() for x in valid_ds], [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))


    # fix encoder parameter
    if fixed_encoder:
        print('fixed parameter of encoder')
        for name, value in net.named_parameters():
            if name.startswith('seq_encoder'):
                value.requires_grad = False

    print("fold", fold_i ,"training data:", len(train_ds), "validated data:", len(valid_ds), 'total data:', len(training_dataset))
    trainloader = data_utils.DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=0, pin_memory=True)
    early_stopping = EarlyStopping(patience=3, verbose=True,
                                   path= saved_model_path + "/fold_" + str(fold_i) + "_best_"+model_name+"_checkpoint.pt")

    
    # get all PSI values from training dataset
    if predict == 'multi' and weigh_samples:
        dw, pos_weight = get_sample_weights(trainloader, device=device)
        print('Dense weight:', dw)
        print('Positive weight:', pos_weight)
    else:
        dw = None
        pos_weight = None

    # Loss functions
    L_expr = nn.MSELoss() #SmoothL1Loss()
    L_splice = nn.BCEWithLogitsLoss()# SmoothL1Loss() # hurdle_loss if predict == 'multi' else nn.SmoothL1Loss()
    learned_loss = True if loss_class is not None else False

    all_params = net.parameters() if not learned_loss else list(net.parameters()) + list(loss_class.parameters())
    optimizer = torch.optim.AdamW(all_params, lr=learning_rate, weight_decay=1e-6)
    net.train()

    # when using the xpu device, optimise the model for the xpu
    if device == 'xpu':
        import intel_extension_for_pytorch as ipex
        print("Using XPU device for training")
        net, optimizer = ipex.optimize(net, optimizer=optimizer)

    print('Model name:', net.name)
    lrs = []
    # last_loss = None
    for epoch in range(epochs):
        net.train()
        print('learning rate:', get_lr(optimizer))
        running_loss = 0
        loss_e = 0
        expression_loss = 0
        splice_loss = 0
        # print('model training mode is:', net.training)
        for data in tqdm(trainloader):
            # print(inputs.size())
            input_pe, input_seg, input_feat, input_dist, y_expr, y_psi, _ = data
            input_pe = input_pe.float().to(device)

            input_seg = input_seg.float().to(device)
            input_feat = input_feat.float().to(device)
            input_dist = input_dist.float().to(device)

            y_expr = y_expr.float().to(device)
            y_psi = y_psi.float().to(device)

            pred_expr, pred_splice_binary, pred_splice, _ = net(input_pe, input_seg, input_feat, input_dist)

            loss_expr = L_expr(pred_expr, y_expr.reshape(pred_expr.shape))

            # if dw is not None:
            #     loss_splice_binary, loss_splice = L_splice(pred_splice_binary, pred_splice, y_psi, loss_type='dense',
            #                                                dw=dw, pos_weight=pos_weight)
            # elif predict == 'multi' and dw is None:
            #     loss_splice_binary, loss_splice = L_splice(pred_splice_binary, pred_splice, y_psi, loss_type='l1')
            # else:
            loss_splice = L_splice(pred_splice, y_psi.reshape(pred_splice.shape))
            # loss_splice_binary = 0.0 # placeholder to keep the code structure consistent

            expression_loss += loss_expr.item()
            splice_loss += loss_splice.item() #if loss_splice is not None else 0.0

            # adding all losses together
            #weight_e, weight_b, weight_s = torch.tensor(1.0, device=device), torch.tensor(0.0, device=device), torch.tensor(0.0, device=device)

            # if learned_loss:
            #     loss, weight_e, weight_b, weight_s = loss_class(expression_loss, loss_splice_binary, splice_loss)

            # elif predict == 'multi':
            #     weight_b, weight_s = torch.tensor(1.0, device=device), torch.tensor(1.0, device=device)
            #
            #     loss = (loss_expr * weight_e) + (loss_splice_binary * weight_b)
            #     if loss_splice is not None:
            #         loss += (loss_splice * weight_s)
            #
            # else:
            #     loss = loss_expr

            # loss_e += loss_expr.item() #+ loss_splice_binary # ((loss_expr.item() * weight_e) + (loss_splice_binary.item() * weight_b))
            # if loss_splice is not None:
            #     loss_e += (loss_splice.item())# * weight_s)

            loss = loss_expr + loss_splice
            running_loss += loss_expr.item() + loss_splice.item() #+ loss_splice_binary.item() if loss_splice_binary is not None else 0.0
            # propagate the loss backward
            loss.backward()

            # update the gradients
            optimizer.step()
            optimizer.zero_grad()

        print(f"[Epoch {epoch + 1}] overall loss: {running_loss/len(trainloader):.5f}, "
              f"expression loss: {expression_loss/len(trainloader):.5f}, " \
              f"splice loss: {splice_loss/len(trainloader):.5f}")
        # print('Training Loss:', loss_e/len(trainloader))
        # print('Loss weights:', weight_e.item(), weight_b.item(), weight_s.item())
        # log_cols = ['Epoch', 'Training_Loss', 'Validation_Loss', 'Validation_PearsonR_allGene',
        #             'Validation_R2_allGene', 'Validation_PearsonR_weGene', 'Validation_R2_weGene', 'Saved?']


        val_mse_all, val_r2_all, val_pr_all, val_loss_total, val_loss_expr, val_loss_splice = validate(net, valid_ds, n_enhancers=n_enhancers, device=device,
                                                       predict=predict)#, loss_weights=(weight_e, weight_b, weight_s))
        val_pr_we, val_r2_we = val_pr_all, val_r2_all
        print('Valdaition R square all:', val_r2_all)
        early_stopping(val_loss_total, net, epoch)
        loss_file.write(f"{fold_i},{epoch+1},{running_loss/len(trainloader)},{expression_loss/len(trainloader)}," \
                        f"{splice_loss/len(trainloader)},{val_mse_all},{val_r2_all},{val_loss_total},{val_loss_expr},{val_loss_splice}\n")
        loss_file.flush()
        if model_logger is not None:
            label_type = net.name.split('.')[-1]
            model_logger.add([fold_i, epoch, running_loss/len(trainloader), val_mse_all, val_pr_all, val_r2_all, 
                              val_pr_we, val_r2_we, early_stopping.counter, label_type])
            # model_logger.save("./EPInfomrer_log/{}.crossValid.log".format(net.name.replace('.'+label_type, '')))
        if early_stopping.early_stop:
            print("Early stopping")
            break
    return lrs

def validate(net, valid_ds,  net_type = 'seq_feat_dist', n_enhancers=50, batch_size=16, device = 'cuda', 
             predict='multi', loss_weights=(1.0, 1.0, 1.0)):
    validloader = data_utils.DataLoader(valid_ds, batch_size=batch_size, pin_memory=True, num_workers=0)
    net.eval()
    L_expr = nn.MSELoss() # SmoothL1Loss()
    L_psi = nn.BCEWithLogitsLoss() # hurdle_loss if predict == 'multi' else nn.SmoothL1Loss()

    with torch.no_grad():
        preds = []
        actual = []
        preds_psi = []
        actual_psi = []
        loss_e = 0
        min_max_expr = [0,0]
        min_max_psi = [0,0]
        expression_loss = 0
        splice_loss = 0
        for data in validloader:
            # print(inputs.size())
            input_pe, input_seg, input_feat, input_dist, y_expr, y_psi, _ = data
            input_pe = input_pe.float().to(device)
            input_seg = input_seg.float().to(device)
            input_feat = input_feat.float().to(device)
            # input_dist = input_dist.long().to(device)
            input_dist = input_dist.float().to(device)
            # print(input_dist.shape, input_dist)
            # input_PEmask = ~(input_PE.sum(-1).sum(-1) > 0).bool().to(device)
            y_expr = y_expr.float().to(device)
            y_psi = y_psi.float().to(device)
            # print(input_P.shape, input_E.shape, input_Emask.shape)
            pred_expr, pred_splice_binary, pred_splice, _ = net(input_pe, input_seg, input_feat, input_dist)

            outputs = list(pred_expr.flatten().cpu().detach().numpy())
            labels = list(y_expr.flatten().cpu().detach().numpy())

            outputs_psi = list(torch.sigmoid(pred_splice).flatten().cpu().detach().numpy())
            labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            # if predict == 'multi':
            #     outputs_psi = combine_hurdle_outputs(pred_splice_binary, pred_splice)
            #     labels_psi = list(y_psi.flatten().cpu().detach().numpy())
            # else:
            #     outputs_psi = list(pred_splice.flatten().cpu().detach().numpy())
            #     labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            loss_expr = L_expr(pred_expr, y_expr.reshape(pred_expr.shape))
            expression_loss += loss_expr.item()
            # loss_expr = anti_bias_loss(loss_expr, pred_expr, alpha=0.3) # anti-bias loss
            # if predict == 'multi':
            #     loss_splice_binary, loss_splice = L_psi(pred_splice_binary, pred_splice, y_psi, loss_type='mse')
            #     loss_e += ((loss_expr.item() * loss_weights[0]) + (loss_splice_binary.item() * loss_weights[1]))
            #     if loss_splice is not None:
            #         loss_e += (loss_splice.item() * loss_weights[2])
            #     else:
            #         # placeholder to keep the code structure consistent, won't be used in the loss calculation
            #         loss_splice = 0.0
            # else:
            loss_splice = L_psi(pred_splice, y_psi.reshape(pred_splice.shape))
            splice_loss += loss_splice.item()
                # loss_e += loss_expr.item()
            loss_e += loss_expr + loss_splice

            preds += outputs
            actual += labels
            preds_psi += outputs_psi
            actual_psi += labels_psi

            if min_max_expr[0] > np.min(outputs) or min_max_expr[1] < np.max(outputs):
                min_max_expr = [np.min(outputs), np.max(outputs)]

            if min_max_psi[0] > np.min(outputs_psi) or min_max_psi[1] < np.max(outputs_psi):
                min_max_psi = [np.min(outputs_psi), np.max(outputs_psi)]

    r2_value = r2_score(actual, preds)
    peasonr, _ = stats.pearsonr(preds, actual)
    mse = mean_squared_error(preds, actual)
    print(f"[Val] overall loss: {loss_e / len(validloader):.5f}, "
          f"expression loss: {expression_loss / len(validloader):.5f}, " \
          f"splice loss: {splice_loss / len(validloader):.5f}")
    print("\n### Validation ### TPM expresion ###")
    # print('### Loss:', loss_expr.item()/len(validloader))#, 'Loss weighted:', loss_expr.item()/len(validloader) * loss_weights[0].item())
    print('### Min-Max Expression:', min_max_expr)
    print("### MSE:", mse, "R²:", r2_value, 'PeasonR:', peasonr)
    # print("###"*20, "\n")

    try:
        r2_value_psi = r2_score(actual_psi, preds_psi)
        peasonr_psi, _ = stats.pearsonr(preds_psi, actual_psi)
        mse_psi = mean_squared_error(preds_psi, actual_psi)
    except ValueError:
        r2_value_psi = 0
        peasonr_psi = 0
        mse_psi = 0
    print("### Validation ### PSI expression ###")
    # print('### Loss:', splice_loss/len(validloader))#, 'Loss weighted:', loss_splice.item()/len(validloader) * loss_weights[2].item())
    print('### Min-Max PSI:', min_max_psi)
    print("### MSE:", mse_psi, "R²:", r2_value_psi, 'PeasonR:', peasonr_psi)
    # print("###"*20)

    # get average r2 
    if predict == 'multi':
        avg_r2 = (r2_value + r2_value_psi) / 2
    else:
        avg_r2 = r2_value

    return mse, avg_r2, peasonr, loss_e / len(validloader), expression_loss / len(validloader), splice_loss / len(validloader)

def test(net, test_ds, fold_i, model_name = None, saved_model_path=None, batch_size=64, device = 'cuda', 
         model_type='best', normals=None, predict='multi'):
    testloader = data_utils.DataLoader(test_ds, batch_size=batch_size, pin_memory=True, num_workers=0)
    # checkpoint = torch.load(saved_model_path + "/fold_" + str(fold_i) + "_"+model_name+"_checkpoint.pt")
    # net.load_state_dict(checkpoint['model_state_dict'])
    # except:
    # net = nn.DataParallel(net, device_ids=[0,1])
    # net.load_state_dict(checkpoint['model_state_dict'])
    # net.load_state_dict(torch.load("./K562_10crx_models/fold_" + str(fold_i) + "_best_"+model_name+"_checkpoint.pt"))
    # print("Load the best model from fold_" + str(fold_i) + "_"+model_type+"_"+model_name+"_checkpoint.pt", )
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
        uncorr_pred_expr = []
        uncorr_pred_splice = []
        uncorr_actual_expr = []
        uncorr_actual_psi = []

        ensid_list = []
        for data in tqdm(testloader):
            input_pe, input_seg, input_feat, input_dist, y_expr, y_psi, eid = data
            input_pe = input_pe.float().to(device)
            input_seg = input_seg.float().to(device)
            input_feat = input_feat.float().to(device)
            # input_dist = input_dist.long().to(device)
            input_dist = input_dist.float().to(device)
            y_expr = y_expr.float().to(device)
            y_psi = y_psi.float().to(device)
            pred_expr, pred_splice_binary, pred_splice, _ = net(input_pe, input_seg, input_feat, input_dist)

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
            labels = list(y_expr.flatten().cpu().detach().numpy())

            if predict == 'multi':
                outputs_psi = combine_hurdle_outputs(pred_splice_binary, pred_splice)
                labels_psi = list(y_psi.flatten().cpu().detach().numpy())
            else:
                outputs_psi = list(pred_splice.flatten().cpu().detach().numpy())
                labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            preds_psi += outputs_psi
            actual_psi += labels_psi

            preds += outputs
            actual += labels
            ensid_list += eid

            if normals:
                uncorr_pred_expr += list(uncorr_pred_ep.flatten().cpu().detach().numpy())
                uncorr_pred_splice += list(uncorr_pred_sp.flatten().cpu().detach().numpy())
                uncorr_actual_expr += list(uncorr_y_expr.flatten().cpu().detach().numpy())
                uncorr_actual_psi += list(uncorr_y_psi.flatten().cpu().detach().numpy())


    sys.stdout.flush()
    df = pd.DataFrame(index=np.array(ensid_list).flatten())
    df['PredExpr'] = preds
    df['ActualExpr'] = actual
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
