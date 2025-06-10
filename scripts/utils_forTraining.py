# -*- coding: utf-8 -*-

import torch
import os
import numpy as np
import pandas as pd
# torch
import torch.nn as nn
# import torch.nn.functional as F
import torch.optim
import torch.utils.data as data_utils
from torch.utils.data import Subset, Dataset
from denseweight import DenseWeight

import sys
import argparse

from scipy import stats
# import sklearn
from sklearn.metrics import mean_squared_error
from tqdm import tqdm
from sklearn.model_selection import train_test_split
# logging
from tqdm import tqdm
# from model.EPInformer import EPInformer_v2, enhancer_predictor_256bp
import h5py
import glob

def get_lr(optimizer):
    for param_group in optimizer.param_groups:
        return param_group['lr']

def anti_mse_loss(pred, target, alpha=0.3):
    mse = (pred - target) ** 2
    confidence_penalty = pred - pred ** 2
    return (mse + alpha * confidence_penalty)

def dense_loss(pred, target, dw): # from https://link.springer.com/article/10.1007/s10994-021-06023-5
    if dw is None:
        return nn.SmoothL1Loss()(pred, target)
    weight = dw(target.float().cpu().detach().numpy()) 
    loss = anti_mse_loss(pred, target).cpu().detach().numpy()
    weighted_loss = weight * loss
    return torch.tensor(weighted_loss.mean(), dtype=torch.float32, device=pred.device)

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

    def start(self):
        """Begin the recording process."""

        self.data = {name: [] for name in self.names}

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

def train(net, training_dataset, fold_i, saved_model_path='../models', learning_rate=1e-4, model_logger=None, fixed_encoder = False, 
          n_enhancers = 50, valid_dataset = None, model_name = '', batch_size = 64, device = 'cuda', stratify=None, 
          class_weight=None, EPOCHS=100, valid_size=1000, predict='multi'):
    if not os.path.exists(saved_model_path):
        os.mkdir(saved_model_path)
    if not os.path.exists(saved_model_path + "/losses.csv"):
        loss_file = open(saved_model_path + "/losses.csv", "w")
        loss_file.write("fold,epoch,training_loss,expresssion_loss,splice_loss,validation_mse,validation_r2\n")
    else:
        loss_file = open(saved_model_path + "/losses.csv", "a")

    if valid_dataset is not None:
        train_ds = training_dataset
        valid_ds = valid_dataset
    else:
        train_idx, val_idx = train_test_split(list(range(len(training_dataset))), test_size=valid_size, 
                                              shuffle=True, random_state=66, stratify=stratify)
        train_ds = Subset(training_dataset, train_idx)
        valid_ds = Subset(training_dataset, val_idx)

    # fix encoder parameter
    if fixed_encoder:
        print('fixed parameter of encoder')
        for name, value in net.named_parameters():
            if name.startswith('seq_encoder'):
                value.requires_grad = False

    print("fold", fold_i ,"training data:", len(train_ds), "validated data:", len(valid_ds), 'total data:', len(training_dataset))
    trainloader = data_utils.DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=5, pin_memory=True)
    early_stopping = EarlyStopping(patience=3,
               verbose=True, path= saved_model_path + "/fold_" + str(fold_i) + "_best_"+model_name+"_checkpoint.pt")

    dw = DenseWeight(alpha=0.6)
    # get all PSI values from training dataset
    if predict == 'multi':
        all_psi = []
        for data in trainloader:
            _, _, _, _, _, y_psi, _ = data
            all_psi.extend(y_psi.flatten().cpu().numpy())
        dw.fit(np.array(all_psi))
    else:
        dw = None
    L_expr = nn.SmoothL1Loss()
    L_splice = nn.SmoothL1Loss()
    loss_weights = (1.0, 1.0) if predict == 'multi' else (1.0, 0.0)
    optimizer = torch.optim.AdamW(net.parameters(), lr=learning_rate, weight_decay=1e-6)
    print('Model name:', net.name)
    lrs = []
    # last_loss = None
    net.train()
    for epoch in range(EPOCHS):
        net.train()
        print('learning rate:', get_lr(optimizer))
        running_loss = 0
        loss_e = 0
        expression_loss = 0
        splice_loss = 0
        # print('model training mode is:', net.training)
        for data in tqdm(trainloader):
            # print(inputs.size())
            optimizer.zero_grad()
            input_PE, input_seg, input_feat, input_dist, y_expr, y_psi, eid = data
            input_PE = input_PE.float().to(device)

            input_seg = input_seg.float().to(device)
            input_feat = input_feat.float().to(device)
            input_dist = input_dist.float().to(device)

            y_expr = y_expr.float().to(device)
            y_psi = y_psi.float().to(device)


            pred_expr, pred_splice, _ = net(input_PE, input_seg, input_feat, input_dist)
            # check devices for all tensors
            loss_expr = L_expr(pred_expr, y_expr)
            if dw is not None:
                loss_splice = dense_loss(pred_splice, y_psi, dw)
            else:
                loss_splice = L_splice(pred_splice, y_psi) # splicing loss here
            loss_e += ((loss_expr.item() * loss_weights[0]) + (loss_splice.item() * loss_weights[1]))
            expression_loss += loss_expr.item()
            splice_loss += loss_splice.item()

            loss = loss_expr# + loss_intensity + loss_contact
            if predict == 'multi':
                loss = (loss_expr * loss_weights[0]) + (loss_splice * loss_weights[1])
            # propagate the loss backward
            loss.backward()
            # update the gradients
            optimizer.step()
            running_loss += loss.item()

        print('[Epoch %d] loss: %.9f' % (epoch + 1, running_loss/len(trainloader)))
        print('Training Loss:', loss_e/len(trainloader))
        # log_cols = ['Epoch', 'Training_Loss', 'Validation_Loss', 'Validation_PearsonR_allGene',
        #             'Validation_R2_allGene', 'Validation_PearsonR_weGene', 'Validation_R2_weGene', 'Saved?']


        val_mse_all, val_r2_all, val_pr_all = validate(net, valid_ds, n_enhancers=n_enhancers, device=device)
        val_r2 = val_r2_all
        val_pr_wE, val_r2_wE = val_pr_all, val_r2_all
        print('Valdaition R square all:', val_r2_all)
        early_stopping(-val_r2, net, epoch)
        loss_file.write(f"{fold_i},{epoch+1},{running_loss/len(trainloader)},{expression_loss/len(trainloader)}," \
                        f"{splice_loss/len(trainloader)},{val_mse_all},{val_r2_all}\n")
        loss_file.flush()
        if model_logger is not None:
            label_type = net.name.split('.')[-1]
            model_logger.add([fold_i, epoch, running_loss/len(trainloader), val_mse_all, val_pr_all, val_r2_all, 
                              val_pr_wE, val_r2_wE, early_stopping.counter, label_type])
            # model_logger.save("./EPInfomrer_log/{}.crossValid.log".format(net.name.replace('.'+label_type, '')))
        if early_stopping.early_stop:
            print("Early stopping")
            break
    return lrs

def validate(net, valid_ds,  net_type = 'seq_feat_dist', n_enhancers=50, batch_size=16, device = 'cuda'):
    validloader = data_utils.DataLoader(valid_ds, batch_size=batch_size, pin_memory=True, num_workers=0)
    net.eval()
    L_expr = nn.MSELoss()
    L_psi = nn.MSELoss()
    
    with torch.no_grad():
        preds = []
        actual = []
        preds_psi = []
        actual_psi = []
        loss_e = 0
        for data in validloader:
            # print(inputs.size())
            input_PE, input_seg, input_feat, input_dist, y_expr, y_psi, eid = data
            input_PE = input_PE.float().to(device)
            input_seg = input_seg.float().to(device)
            input_feat = input_feat.float().to(device)
            # input_dist = input_dist.long().to(device)
            input_dist = input_dist.float().to(device)
            # print(input_dist.shape, input_dist)
            # input_PEmask = ~(input_PE.sum(-1).sum(-1) > 0).bool().to(device)
            y_expr = y_expr.float().to(device)
            y_psi = y_psi.float().to(device)
            # print(input_P.shape, input_E.shape, input_Emask.shape)
            pred_expr, pred_splice, _ = net(input_PE, input_seg, input_feat, input_dist)

            outputs = list(pred_expr.flatten().cpu().detach().numpy())
            labels = list(y_expr.flatten().cpu().detach().numpy())

            outputs_psi = list(pred_splice.flatten().cpu().detach().numpy())
            labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            loss_expr = L_expr(pred_expr, y_expr)
            loss_splice = L_psi(pred_splice, y_psi)
            loss_e += loss_expr.item() + loss_splice.item()

            preds += outputs
            actual += labels
            preds_psi += outputs_psi
            actual_psi += labels_psi

    slope, intercept, r_value, p_value, std_err = stats.linregress(preds, actual)
    peasonr, pvalue = stats.pearsonr(preds, actual)
    mse = mean_squared_error(preds, actual)
    print("### Validation ### TPM expresion ###")
    print('### Loss:', loss_expr.item()/len(validloader))
    print("### MSE:", mse, "R²:", r_value**2, 'PeasonR:', peasonr)
    print("###"*20, "\n")

    try:
        slope, intercept, r_value_psi, p_value, std_err = stats.linregress(preds_psi, actual_psi)
        peasonr_psi, pvalue = stats.pearsonr(preds_psi, actual_psi)
        mse_psi = mean_squared_error(preds_psi, actual_psi)
    except ValueError:
        r_value_psi = 0
        peasonr_psi = 0
        mse_psi = 0
    print("### Validation ### PSI expression ###")
    print('### Loss:', loss_splice.item()/len(validloader))
    print("### MSE:", mse_psi, "R²:", r_value_psi**2, 'PeasonR:', peasonr_psi)
    print("###"*20)

    return mse, r_value**2, peasonr

def test(net, test_ds, fold_i, model_name = None, saved_model_path=None, batch_size=64, device = 'cuda', model_type='best'):
    testloader = data_utils.DataLoader(test_ds, batch_size=batch_size, pin_memory=True, num_workers=0)
    # checkpoint = torch.load(saved_model_path + "/fold_" + str(fold_i) + "_"+model_name+"_checkpoint.pt")
    # net.load_state_dict(checkpoint['model_state_dict'])
    # except:
    # net = nn.DataParallel(net, device_ids=[0,1])
    # net.load_state_dict(checkpoint['model_state_dict'])
    # net.load_state_dict(torch.load("./K562_10crx_models/fold_" + str(fold_i) + "_best_"+model_name+"_checkpoint.pt"))
    # print("Load the best model from fold_" + str(fold_i) + "_"+model_type+"_"+model_name+"_checkpoint.pt", )
    if saved_model_path is not None:
        checkpoint = torch.load(saved_model_path + "/fold_" + str(fold_i) + "_best_"+model_name+"_checkpoint.pt", weights_only=False)
        net.load_state_dict(checkpoint['model_state_dict'])
        print(model_name,'loaded!')
        
    net.eval()
    with torch.no_grad():
        preds = []
        actual = []
        preds_psi = []
        actual_psi = []

        ensid_list = []
        for data in tqdm(testloader):
            input_PE, input_seg, input_feat, input_dist, y_expr, y_psi, eid = data
            input_PE = input_PE.float().to(device)
            input_seg = input_seg.float().to(device)
            input_feat = input_feat.float().to(device)
            # input_dist = input_dist.long().to(device)
            input_dist = input_dist.float().to(device)
            y_expr = y_expr.float().to(device)
            y_psi = y_psi.float().to(device)
            pred_expr, pred_splice, _ = net(input_PE, input_seg, input_feat, input_dist)

            outputs = list(pred_expr.flatten().cpu().detach().numpy())
            labels = list(y_expr.flatten().cpu().detach().numpy())

            outputs_psi = list(pred_splice.flatten().cpu().detach().numpy())
            labels_psi = list(y_psi.flatten().cpu().detach().numpy())

            preds_psi += outputs_psi
            actual_psi += labels_psi

            preds += outputs
            actual += labels
            ensid_list += eid


    slope, intercept, r_value, p_value, std_err = stats.linregress(preds, actual)
    peasonr, pvalue = stats.pearsonr(preds, actual)
    mse = mean_squared_error(preds, actual)
    # print(fold %s test sequence: %0.3f' % (fold_i, r_value**2))
    sys.stdout.flush()
    df = pd.DataFrame(index=np.array(ensid_list).flatten())
    df['PredExpr'] = preds
    df['ActualExpr'] = actual
    if len(preds) == len(preds_psi):
        df['PredPsi'] = preds_psi
        df['ActualPsi'] = actual_psi
    df['fold_idx'] = fold_i

    if saved_model_path is not None:
        df.to_csv(saved_model_path + "/fold_" + str(fold_i) + "_"+ model_name + "_predictions.csv")
    return df
