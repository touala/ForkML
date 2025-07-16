#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from tqdm import tqdm

# For 3d plot
# import json
import glob
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import hpbandster.core.result as hpres

from tensorflow.keras.optimizers import SGD, Adam
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from unet import custom_unet
from metrics import mean_iou, jaccard_distance, MeanX0Distance, compute_x0_distance
from tensorflow.keras.metrics import MeanIoU

#################### Functions ####################

# Currently unused
def get_callbacks():
    early_stopping = EarlyStopping(
        monitor='val_mean_iou', 
        patience=20, 
        verbose=1, 
        mode='max', 
        restore_best_weights=True
    )
    reduce_lr = ReduceLROnPlateau(
        monitor='val_loss', 
        factor=0.2, 
        patience=5, 
        verbose=1, 
        mode='min', 
        min_delta=0.0001, 
        cooldown=0, 
        min_lr=0
    )
    return [early_stopping, reduce_lr]

def save_history(history, filename):
    # convert the history.history dict to a pandas DataFrame
    hist_df = pd.DataFrame(history.history) 
    # and save to csv:
    with open(filename, mode='w') as f:
        hist_df.to_csv(f)

def rle_with_indices(arr):
    """
    Optimized RLE with NumPy to quickly find runs of identical values.
    """
    if len(arr) == 0:
        return []
    
    # Find where the value changes
    changes = np.where(np.diff(arr) != 0)[0]
    start_indices = np.concatenate(([0], changes + 1))
    end_indices = np.concatenate((changes, [len(arr) - 1]))
    
    # Return list of (value, start_index, end_index)
    rle_list = list(zip(arr[start_indices], start_indices, end_indices))
    
    return rle_list

def update_color_spans(ax, class_idx, class_colors, ymin, ymax, masked_indices=None):
    if masked_indices is not None:
        class_idx[masked_indices] = 3  # Assign the masked value (3) to masked_indices
    
    # Use RLE to get start and end indices of runs
    rle_result = rle_with_indices(class_idx)
    
    for value, start_idx, end_idx in rle_result:
        ax.axvspan(start_idx, end_idx + 1, ymin, ymax, facecolor=class_colors[value], alpha=0.3)

def plot_signals_with_colored_labels(signals, labels, output_filename, predicted_labels=None, with_mask=None, path_valid_id=None):
    if path_valid_id is None:
        path_valid_id = 'valid_read_ids.txt'
    
    with open(path_valid_id, 'r') as file:
        # Read all lines from the file and strip newline characters
        read_ids = [line.strip() for line in file]
    
    if with_mask is not None:
        if with_mask == 'rough':
            masked_value = 0
        elif with_mask == 'clean':
            masked_value = -1
    
    with PdfPages(output_filename) as pdf:
        num_samples = len(signals) # 50
        class_colors = {0: 'red', 1: 'blue', 2: 'beige', 3: 'lightgrey'}
        class_names = {0: 'right', 1: 'left', 2: 'none', 3: 'masked'}
        
        # Using tqdm to add a progress bar to the loop
        for i in tqdm(range(num_samples), desc="Processing signals"):
            # print(f"Processing Sample {i+1}")
            class_idx = np.argmax(labels[i], axis=1)
            contains_class_0_or_1 = (0 in class_idx) or (1 in class_idx)
            if not contains_class_0_or_1 and predicted_labels is not None:
                class_idx_pred = np.argmax(predicted_labels[i], axis=1)
                contains_class_0_or_1 = contains_class_0_or_1 or (0 in class_idx_pred) or (1 in class_idx_pred)
            
            if not contains_class_0_or_1:
                continue  # Skip this signal if neither labels nor predicted_labels contain class 0 or 1
            
            signal_length = len(signals[i])
            fig, ax = plt.subplots(figsize=(10, 3))
            ax.plot(signals[i], label=f"Signal {i + 1}", color="black")
            
            masked_indices = None
            if with_mask is not None:
                masked_indices = np.where(signals[i] == masked_value)[0]  # indices where signal is -1
            
            ymin_labels = 0.5 if predicted_labels is not None else 0
            update_color_spans(ax, class_idx, class_colors, ymin_labels, 1, masked_indices)
            
            if predicted_labels is not None:
                class_idx_pred = np.argmax(predicted_labels[i], axis=1)
                update_color_spans(ax, class_idx_pred, class_colors, 0, 0.5, masked_indices)
                ax.text(0.95, 0.25, 'Predicted', transform=ax.transAxes, ha='center', va='center', fontsize=10)
                
                # Annotate the plot with distances and positions
                distances, positions = compute_x0_distance(labels[i], predicted_labels[i])
                # print(distance    s, positions)
                for j, pos in enumerate(positions):
                    ax.annotate(f"{distances[j]}", (pos, 0.9), textcoords="offset points", xytext=(0, 0), ha='center', color='green', fontsize=8)
            
            ax.text(0.95, 0.75, 'Annotated', transform=ax.transAxes, ha='center', va='center', fontsize=10)
            ax.set_title(f"Sample {i + 1} - {read_ids[i]}")
            ax.set_ylim([0, 1])
            ax.set_xlim([0, signal_length])
            legend_handles = [Patch(facecolor=class_colors[key], label=class_names[key], alpha=0.3) for key in class_colors]
            ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
            
            plt.subplots_adjust(left=0.01, right=0.95, top=0.85, bottom=0.2)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

# https://keras.io/guides/keras_tuner/getting_started/
def call_existing_code(nb_filters, nb_layers, sz_kernel, masking_type, lr, eps, use_batch_norm, dropout, dropout_change_per_layer):
    model = custom_unet(
        input_shape=(None, 1),
        num_classes=3,
        activation="relu",
        use_batch_norm=use_batch_norm,
        dropout=dropout,
        dropout_change_per_layer=dropout_change_per_layer,
        dropout_type="spatial",
        use_dropout_on_upsampling=False,
        filters=nb_filters,
        num_layers=nb_layers,
        kernel_size=sz_kernel,
        output_activation="softmax",
        include_masking=masking_type
    )
    
    opt = Adam(learning_rate=lr, epsilon=eps)
    model.compile(
        optimizer = opt,
        loss = "categorical_crossentropy", # alt tf.keras.losses.CategoricalCrossentropy
        metrics = [mean_iou, MeanIoU(num_classes=3), MeanX0Distance()] # jaccard_distance, MeanIoU(num_classes=3, ignore_class=2) not working
    )
    return(model)

def build_model(hp):
    nb_filters = hp.Int("nb_filters", min_value=4, max_value=64, step=1)
    nb_layers = hp.Int("nb_layers", min_value=4, max_value=7, step=1) # OOM at 10
    sz_kernel = hp.Int("sz_kernel", min_value=6, max_value=32, step=1)
    # activation = hp.Choice("activation", ["relu", "tanh"])
    # dropout = hp.Boolean("dropout")
    lr = hp.Float("lr", min_value=1e-4, max_value=1e-2, sampling="log")
    eps = hp.Float("eps", min_value=1e-7, max_value=1e-4, sampling="log")
    # call existing model-building code with the hyperparameter values.
    model = call_existing_code(
        nb_filters=nb_filters, nb_layers=nb_layers, sz_kernel=sz_kernel, lr=lr, eps=eps #  activation=activation, dropout=dropout
    )
    return model

def hbtuning(hyperband_iterations, max_epochs, factor):
    import math
    return round(hyperband_iterations*max_epochs*pow(math.log(max_epochs,factor),2))

def summarize_bohb_results(path_dir):
    result = hpres.logged_results_to_HBS_result(path_dir)
    
    # get all executed runs
    all_runs = result.get_all_runs()
    
    # get the 'dict' that translates config ids to the actual configurations
    id2conf = result.get_id2config_mapping()
    
    data = []
    for run in all_runs:  # Assuming all_runs is the list you provided
        # Extracting details from each run
        if run.info is None:
            entry = {
                'config_id': run.config_id,
                'budget': run.budget,
                'best_epoch': None,
                'loss': run.loss,
                # 'submitted': run.time_stamps['submitted'],
                'start': run.time_stamps['started'],
                'end': run.time_stamps['finished'],
                'train_loss': None,
                'train_acc': None,
                'train_acc2': None,
                'train_dist': None,
                'val_loss': None,
                'val_acc': None,
                'val_acc2': None,
                'val_dist': None,
                'nb_parameters': None
            }
        else:
            entry = {
                'config_id': run.config_id,
                'budget': run.budget,
                'best_epoch': run.info['best epoch'],
                'loss': run.loss,
                # 'submitted': run.time_stamps['submitted'],
                'start': run.time_stamps['started'],
                'end': run.time_stamps['finished'],
                'train_loss': run.info['train loss'],
                'train_acc': run.info['train accuracy'],
                'train_acc2': run.info['train accuracy2'],
                'train_dist': run.info['train x0 dist'],
                'val_loss': run.info['validation loss'],
                'val_acc': run.info['validation accuracy'],
                'val_acc2': run.info['validation accuracy2'],
                'val_dist': run.info['validation x0 dist'],
                'nb_parameters': run.info['number of parameters']
            }
        data.append(entry)
    
    configs = []
    for key, value in id2conf.items():
        flattened_dict = {'config_id': key}
        flattened_dict.update(value['config'])
        flattened_dict.update(value['config_info'])
        configs.append(flattened_dict)
    
    hpt_results = pd.merge(pd.DataFrame(configs), pd.DataFrame(data), on='config_id', how='inner')
    
    return hpt_results
