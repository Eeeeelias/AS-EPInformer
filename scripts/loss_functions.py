import torch
import torch.nn as nn


def anti_mse_loss(pred, target, alpha=0.3):
    mse = (pred - target) ** 2
    confidence_penalty = pred - pred ** 2
    return mse + alpha * confidence_penalty


def dense_loss(pred, target, dw): # from https://link.springer.com/article/10.1007/s10994-021-06023-5
    if dw is None:
        return nn.SmoothL1Loss()(pred, target)
    weight = dw(target.float().cpu().detach().numpy()) 
    loss = anti_mse_loss(pred, target, alpha=1).cpu().detach().numpy()
    weighted_loss = weight * loss
    return torch.tensor(weighted_loss.mean(), dtype=torch.float32, device=pred.device)


def hurdle_loss(binary_logits, splicing_pred, target, loss_type='l1', dw=None, pos_weight=None):
    binary_target = (target != 1).float()

    binary_target = binary_target.unsqueeze(-1)
    bce_loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    bce_loss = bce_loss_fn(binary_logits, binary_target)

    # regression part only for values != 1
    regression_mask = (target != 1)
    if regression_mask.any():
        if loss_type == 'l1':
            regression_loss_fn = nn.SmoothL1Loss()
        elif loss_type == 'mse':
            regression_loss_fn = nn.MSELoss()
        elif loss_type == 'dense': 
            regression_loss_fn = lambda x, y: dense_loss(x, y, dw)
        else:
            raise ValueError(f"Unsupported loss type: {loss_type}")

        regression_loss = regression_loss_fn(splicing_pred[regression_mask], target[regression_mask])
    else:
        regression_loss = None

    return bce_loss, regression_loss


def anti_bias_loss(mse_loss, pred, alpha=0.3):
    variance_penalty = -torch.var(pred)
    return mse_loss + alpha * variance_penalty
