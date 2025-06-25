import math
import torch
import csv
import numpy as np
import matplotlib.pyplot as plt

from EPInformer.models_multi import EPInformer_v2, enhancer_predictor_256bp, WeightedLoss

# device = torch.device('xpu' if torch.xpu.is_available() else 'cpu')
device = torch.device('cpu')


pretrained_convNet = enhancer_predictor_256bp()
pt_model_name = 'K562_seq2activityLog2_leaveChrOut_combinedRS_2bins_bs64_H3K27ac_adamW_erisxdl_r0'
checkpoint = torch.load(f"./trained_models/pretrained_enhancer_encoder/fold_1_best_{pt_model_name}_checkpoint.pt", 
                                weights_only=False, map_location=device)

model = EPInformer_v2(n_encoder=3, pre_trained_encoder=pretrained_convNet.encoder,
                              n_enhancer=60, out_dim=64, n_extraFeat=3, device=device, 
                              exon_data=True, head=8).to(device)

checkpoint_path= "trained_models/2025-06-19-10/fold_1_best_EPInformer-PE-Activity-HiC.preTrainedConv.4base.64dim." \
                 "3Trans.8head.TrueBN.TrueLN.TrueFeat.3extraFeat.60enh.Trueexon.K562.RNA_checkpoint.pt"
model_checkpoint = torch.load(checkpoint_path, weights_only=False, map_location=device)

# print(model_checkpoint_alt)
model_name = checkpoint_path.split('/')[1]
model.load_state_dict(model_checkpoint['model_state_dict'])


model.eval()

for name, param in model.named_parameters():
    print(f"Parameter: {name}, shape: {param.shape}")


def describe_weights(tensor, name, as_csv=False):
    with torch.no_grad():
        if as_csv:
            return {'name': name, 'shape': tensor.shape, 'mean': tensor.mean().item(), 'std': tensor.std().item(),
             'min': tensor.min().item(), 'max': tensor.max().item(), 'norm': tensor.norm().item(),
             'sparsity': (tensor == 0).float().mean().item() * 100}

        print(f"=== {name} ===")
        print(f"Shape: {tensor.shape}")
        print(f"Mean: {tensor.mean().item():.4f}")
        print(f"Std: {tensor.std().item():.4f}")
        print(f"Min: {tensor.min().item():.4f}")
        print(f"Max: {tensor.max().item():.4f}")
        print(f"Norm: {tensor.norm().item():.4f}")
        print(f"Sparsity (% zeros): {(tensor == 0).float().mean().item() * 100:.2f}%")
        print()


def plot_weight_distributions(model, filter_str="pTo", xlim=None, auto_scale=False, percentile=99):
    # Collect matching parameters
    filtered_params = [(name, param) for name, param in model.named_parameters() if filter_str in name]
    n = len(filtered_params)

    if n == 0:
        print(f"No parameters found containing '{filter_str}'")
        return
    
    if auto_scale:
        all_weights = torch.cat([p.data.flatten().cpu() for _, p in filtered_params])
        limit = np.percentile(np.abs(all_weights.numpy()), percentile)
        xlim = (-limit, limit)

    # Define grid size (auto square-ish)
    cols = min(4, n)
    rows = math.ceil(n / cols)

    # Create subplots
    fig, axs = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3))
    axs = axs.flatten() if n > 1 else [axs]

    for i, (name, param) in enumerate(filtered_params):
        tensor = param.data.flatten().cpu().numpy()
        axs[i].hist(tensor, bins=100, color='steelblue', alpha=0.7)
        axs[i].set_title(name, fontsize=9)
        axs[i].set_xlabel('Weight')
        axs[i].set_ylabel('Freq')
        if xlim:
            axs[i].set_xlim(xlim)
    
    # Turn off empty subplots
    for j in range(i + 1, len(axs)):
        axs[j].axis('off')

    plt.tight_layout()
    plt.savefig(f"images/{model_name}_weight_distributions_{filter_str}.png", dpi=150)
    plt.show()

plot_weight_distributions(model, filter_str="attn", auto_scale=True)


# Save weights summary to CSV
csv_writer = csv.writer(open(f"data/{model_name}_weights_summary.csv", "w", newline='\n'))
csv_writer.writerow(['name', 'shape', 'mean', 'std', 'min', 'max', 'norm', 'sparsity'])
for name, param in model.named_parameters(): 
    out = describe_weights(param.data, name, as_csv=True)
    if out is None:
        print(f"Skipping {name} as it is None")
        continue
    csv_writer.writerow([out['name'], out['shape'], out['mean'], out['std'], out['min'],
                            out['max'], out['norm'], out['sparsity']])
