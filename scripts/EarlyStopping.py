import numpy as np
import torch

class EarlyStoppingMulti:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, expr_patience=2, splice_patience=2, verbose=False, delta=0, path=None):
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
        self.expr_patience = expr_patience
        self.splice_patience = splice_patience
        self.verbose = verbose
        self.expr_counter = 0
        self.splice_counter = 0
        self.best_expr_score = None
        self.best_splice_score = None
        self.expr_early_stop = False
        self.splice_early_stop = False
        self.val_loss_min = np.inf
        self.val_loss_expr = np.inf
        self.val_loss_splice = np.inf
        self.delta = delta
        self.path = path

    def __call__(self, val_loss_expr, val_loss_splice, model, epoch_i):
        expr_score = -val_loss_expr  # We want to minimize the expression loss
        splice_score = -val_loss_splice  # We want to minimize the splice loss

        if self.best_expr_score is None or self.best_splice_score is None:
            self.best_expr_score = expr_score
            self.best_splice_score = splice_score
            self.save_checkpoint(val_loss_expr, val_loss_splice, model, epoch_i)

        if expr_score < self.best_expr_score + self.delta and not self.expr_early_stop:
            self.expr_counter += 1
            print(f"EarlyStopping expr_counter: {self.expr_counter} out of {self.expr_patience}, " \
                    f"expr_score: {self.best_expr_score:.4f}, splice_score: {self.best_splice_score:.4f}")
            if self.expr_counter >= self.expr_patience:
                self.expr_early_stop = True
        
        if splice_score < self.best_splice_score + self.delta and not self.splice_early_stop:
            self.splice_counter += 1
            print(f"EarlyStopping splice_counter: {self.splice_counter} out of {self.splice_patience}, " \
                  f"expr_score: {self.best_expr_score:.4f}, splice_score: {self.best_splice_score:.4f}")
            if self.splice_counter >= self.splice_patience:
                self.splice_early_stop = True

        if expr_score > self.best_expr_score + self.delta or splice_score > self.best_splice_score + self.delta:
            if expr_score > self.best_expr_score + self.delta and not self.expr_early_stop:
                self.best_expr_score = expr_score
                self.expr_counter = 0
            if splice_score > self.best_splice_score + self.delta and not self.splice_early_stop:
                self.best_splice_score = splice_score
                self.splice_counter = 0

            self.save_checkpoint(val_loss_expr, val_loss_splice, model, epoch_i)



    def save_checkpoint(self, val_loss_expr, val_loss_splice, model, epoch_i):
        '''Saves model when validation loss decrease.'''
        val_loss = val_loss_expr + val_loss_splice
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f},' \
                  f'expr: {self.val_loss_expr:.6f} --> {val_loss_expr:.6f},' \
                  f'splice: {self.val_loss_splice:.6f} --> {val_loss_splice:.6f}).  Saving model ...')

        if self.path is not None:
            torch.save({
                        'epoch': epoch_i,
                        'model_state_dict': model.state_dict(),
                        'loss': val_loss,
                        'val_loss_expr': val_loss_expr,
                        'val_loss_splice': val_loss_splice,
                        },
                        self.path)
            print('Saving ckpt at', self.path)
        # torch.save(model.state_dict(), self.path)
        self.val_loss_min = val_loss
        self.val_loss_expr = val_loss_expr
        self.val_loss_splice = val_loss_splice


class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, patience=3, verbose=False, delta=0, path=None):
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
        score = -val_loss.cpu().detach().numpy()  # We want to minimize the validation loss
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

        if self.path is not None:
            torch.save({
                        'epoch': epoch_i,
                        'model_state_dict': model.state_dict(),
                        'loss': val_loss,
                        },
                        self.path)
            print('Saving ckpt at', self.path)
        # torch.save(model.state_dict(), self.path)
        self.val_loss_min = val_loss
