#!/usr/bin/env python3

print("Loading dependencies")

#################### Dependencies ####################
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from tqdm import tqdm

#################### Main ####################
print("Read arguments")
file_train_sigs_bal0 = sys.argv[1]
file_train_labs_bal0 = sys.argv[2]
file_train_sigs_bal3 = sys.argv[3]
file_train_labs_bal3 = sys.argv[4]
file_train_sigs_nobal = sys.argv[5]
file_train_labs_nobal = sys.argv[6]
file_valid_sigs = sys.argv[7]
file_valid_labs = sys.argv[8]
masking_type = sys.argv[9]
with_hp_tuning = sys.argv[10].lower() == "true"
out_dir = sys.argv[11]
tmp_dir = sys.argv[12]
nb_threads = sys.argv[13]
bohb_verbose = 0

os.environ["OMP_NUM_THREADS"] = f'{nb_threads}'

print("Load tensorflow")
import tensorflow as tf
from tensorflow.keras.optimizers import SGD, Adam
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from unet import custom_unet
from training_function import save_history, rle_with_indices, update_color_spans, plot_signals_with_colored_labels, call_existing_code, build_model, hbtuning, summarize_bohb_results
from metrics import mean_iou, jaccard_distance, MeanX0Distance, compute_x0_distance
from tensorflow.keras.metrics import MeanIoU
from tensorflow.keras.models import load_model

import keras_tuner

# Load data
print("Load data")
train_sigs_bal0 = np.load(file_train_sigs_bal0)
train_labs_bal0 = np.load(file_train_labs_bal0)
train_sigs_bal3 = np.load(file_train_sigs_bal3)
train_labs_bal3 = np.load(file_train_labs_bal3)
train_sigs_nobal = np.load(file_train_sigs_nobal)
train_labs_nobal = np.load(file_train_labs_nobal)
valid_sigs = np.load(file_valid_sigs)
valid_labs = np.load(file_valid_labs)

print(f'Size datasets:\n\t(Train fork only: {len(train_sigs_bal0)})\n\tTrain mostly fork (at least 30%): {len(train_sigs_bal3)}\n\t(Train natural content: {len(train_sigs_nobal)})\n\tValidation: {len(valid_sigs)}\n')

if(with_hp_tuning):
    print("Hyperparameters optimization")
    import pickle
    from training_function import save_history, rle_with_indices, update_color_spans, plot_signals_with_colored_labels, call_existing_code, build_model, hbtuning, summarize_bohb_results
    
    from keras.utils import plot_model
    
    from hpbandster.core.worker import Worker
    from hpbandster.optimizers import BOHB
    import hpbandster.core.nameserver as hpns
    import hpbandster.core.result as hpres
    import hpbandster.visualization as hpvis
    
    from ConfigSpace import ConfigurationSpace, UniformFloatHyperparameter, UniformIntegerHyperparameter, CategoricalHyperparameter
    from ConfigSpace.read_and_write import json as cs_json
    
    import shutil
    
    class KerasWorker(Worker):
        def __init__(self, training_dir, bohb_verbose, tmpdir, **kwargs):
            super().__init__(**kwargs)
            self.training_dir = training_dir
            self.bohb_verbose = bohb_verbose
            self.tmpdir = tmpdir
        
        def compute(self, config, budget, **kwargs):
            model = call_existing_code(
                nb_filters=config['nb_filters'], nb_layers=config['nb_layers'], sz_kernel=config['sz_kernel'], masking_type=masking_type, lr=config['lr'], eps=config['eps'], use_batch_norm=config['use_batch_norm'], dropout=config['dropout'], dropout_change_per_layer=config['dropout_change_per_layer'] #  activation=activation
            )
            config_str = "_".join(f"{k}-{v}" for k, v in config.items())
            
            early_stopping = EarlyStopping(monitor='val_mean_iou', patience=12, verbose=self.bohb_verbose, mode='max') # restore_best_weights not used here
            reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=3, verbose=self.bohb_verbose, mode='min', min_delta=0.0001, cooldown=0, min_lr=0)
            
            if int(budget) == pow(my_eta,4):
                saving_best = True
                saved_model = "best_" + config_str + ".h5"
                checkpoint = ModelCheckpoint(filepath=self.tmpdir + '/' + saved_model, monitor="val_mean_iou", mode="max", save_best_only=True, verbose=self.bohb_verbose) # Save best only
                main_callbacks = [early_stopping, reduce_lr, checkpoint]
            else:
                saving_best = False
                main_callbacks = [early_stopping, reduce_lr]
            
            # Use budget as the number of epochs
            train_history = model.fit(train_sigs_bal0, train_labs_bal0, epochs=int(budget), verbose=self.bohb_verbose, validation_data=(valid_sigs, valid_labs), callbacks=main_callbacks) # batch_size=
            save_history(train_history, path_output + "/history_" + config_str + "_" + str(budget) + ".csv")
            idx_best = int(np.argmax(train_history.history['val_mean_iou']))
            
            if saving_best: # Retrieve model from tmpdir
                shutil.move(os.path.join(self.tmpdir, saved_model), os.path.join(self.training_dir, saved_model))
                        
            return ({
                'loss': 1 - float(np.max(train_history.history['val_mean_iou'])), # iou? remember: HpBandSter always minimizes!
                'info': {
                    'state': 'SUCCESS',
                    'train loss': train_history.history['loss'][idx_best],
                    'train accuracy': train_history.history['mean_iou'][idx_best],
                    'train accuracy2': train_history.history['mean_io_u'][idx_best],
                    'train x0 dist': train_history.history['mean_x0_distance'][idx_best],
                    'validation loss': train_history.history['val_loss'][idx_best],
                    'validation accuracy': train_history.history['val_mean_iou'][idx_best],
                    'validation accuracy2': train_history.history['val_mean_io_u'][idx_best],
                    'validation x0 dist': train_history.history['val_mean_x0_distance'][idx_best],
                    'max epoch': len(train_history.history['loss']),
                    'best epoch': int(idx_best + 1), # convert row_idx to epoch
                    'best accuracy': float(np.max(train_history.history['val_mean_iou'])),
                    'number of parameters': model.count_params()
                }
            })
        @staticmethod
        def get_configspace():
            config_space = ConfigurationSpace()
            config_space.add_hyperparameter(UniformIntegerHyperparameter('nb_filters', lower=4, upper=128))
            config_space.add_hyperparameter(UniformIntegerHyperparameter('nb_layers', lower=4, upper=6))
            config_space.add_hyperparameter(UniformIntegerHyperparameter('sz_kernel', lower=6, upper=32))
            config_space.add_hyperparameter(UniformFloatHyperparameter('lr', lower=1e-5, upper=1e-2, log=True))
            config_space.add_hyperparameter(UniformFloatHyperparameter('eps', lower=1e-7, upper=1e-4, log=True))
            config_space.add_hyperparameter(CategoricalHyperparameter('use_batch_norm', choices=[True, False]))
            config_space.add_hyperparameter(UniformFloatHyperparameter('dropout', lower=0, upper=0.9))
            config_space.add_hyperparameter(UniformFloatHyperparameter('dropout_change_per_layer', lower=0, upper=0.1))
            
            configspace_json = cs_json.write(config_space)
            with open(path_output + "/configspace.json", 'w') as f:
                print(path_output)
                f.write(configspace_json)
            
            return(config_space)
    
    tmpdir = tmp_dir
    path_output = f"{out_dir}BOHB_initial"
    path_output_new = path_output
    cpt = 1
    while os.path.exists(path_output_new):
        path_output_new = path_output + "_" + str(cpt)
        cpt += 1
    
    result_logger = hpres.json_result_logger(directory=path_output_new, overwrite=False)
    
    previous_result = None
    if cpt > 1:
        try:
            previous_result = hpres.logged_results_to_HBS_result(path_output)
            path_output = path_output_new
        except Exception as e:
            print(f"Failed to load previous results: {e}")
            previous_result = None
    
    NS = hpns.NameServer(run_id='example', host='localhost', port=None, working_directory=path_output)
    NS.start()
    
    worker = KerasWorker(training_dir=path_output, tmpdir=tmpdir, nameserver='localhost', run_id='example', bohb_verbose=bohb_verbose)
    worker.run(background=True)
    
    my_eta = 3 # 3
    bohb = BOHB( # https://blog.dataiku.com/a-slightly-better-budget-allocation-for-hyperband
        configspace=worker.get_configspace(),
        run_id='example',
        nameserver='localhost',
        min_budget=pow(my_eta,1),
        max_budget=pow(my_eta,4),
        eta=my_eta,
        result_logger=result_logger,
        min_points_in_model = 30,
        previous_result = previous_result
    )
    res = bohb.run(n_iterations=20) # at 20 eta=3 + max 81 -> 230 configs == 324 trainings
    
    bohb.shutdown(shutdown_workers=True)
    NS.shutdown()
    
    hpt_results = summarize_bohb_results(path_output)
    max_config = hpt_results.loc[hpt_results.dropna().groupby('config_id')['val_acc'].idxmax()].sort_values('val_acc').tail(1)
    
    selected_config = {
        'dropout' : max_config['dropout'].item(),
        'dropout_change_per_layer' : max_config['dropout_change_per_layer'].item(),
        'eps' : max_config['eps'].item(),
        'lr' : max_config['lr'].item(),
        'nb_filters' : max_config['nb_filters'].item(),
        'nb_layers' : max_config['nb_layers'].item(),
        'sz_kernel' : max_config['sz_kernel'].item(),
        'use_batch_norm' : max_config['use_batch_norm'].item(),
        'masking_type' : masking_type
    }
    print(f"Training with max. config")
    pd.set_option('display.max_columns', None)
    print(selected_config)
else:
    path_output = f"{out_dir}default_initial"
    
    selected_config = {
        'dropout' : 0.261876500579688,
        'dropout_change_per_layer' : 0.07012775147234629,
        'eps' : 2.28787402711115e-07,
        'lr' : 0.0007597076551517235,
        'nb_filters' : 114,
        'nb_layers' : 4,
        'sz_kernel' : 23,
        'use_batch_norm' : True,
        'masking_type' : "clean"
    }
    print(f"Training with default config")
    pd.set_option('display.max_columns', None)
    print(selected_config)


dropout = selected_config['dropout']
dropout_change_per_layer = selected_config['dropout_change_per_layer']
eps = selected_config['eps']
lr = selected_config['lr']
nb_filters = selected_config['nb_filters']
nb_layers = selected_config['nb_layers']
sz_kernel = selected_config['sz_kernel']
use_batch_norm = selected_config['use_batch_norm']
masking_type = selected_config['masking_type']

model_file = os.path.join(path_output, f"model_final.h5")
history_file = os.path.join(path_output, f"history_final.json")

# Build the model
neuralnet = custom_unet(
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
    include_masking=masking_type  # Set masking type if needed
)

# Compile the model
opt = Adam(learning_rate=lr, epsilon=eps)
neuralnet.compile(
    optimizer=opt,
    loss="categorical_crossentropy",
    metrics=[mean_iou, MeanIoU(num_classes=3), MeanX0Distance()]
)

# Callbacks
early_stopping = EarlyStopping(monitor='val_mean_iou', patience=12, verbose=1, mode='max', restore_best_weights=True)
reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=3, verbose=1, mode='min', min_delta=0.0001, cooldown=0, min_lr=0)
checkpoint = ModelCheckpoint(filepath=model_file, monitor='val_mean_iou', verbose=1, mode='max', save_best_only=True)
main_callbacks = [early_stopping, reduce_lr, checkpoint]

history = neuralnet.fit(train_sigs_bal3, train_labs_bal3, epochs=60, validation_data=(valid_sigs, valid_labs), callbacks=main_callbacks)
save_history(history, history_file)

# Evaluate the model after training
val_score = neuralnet.evaluate(valid_sigs, valid_labs)
print(f"Validation score for max. config: {val_score}")

exit(0)
