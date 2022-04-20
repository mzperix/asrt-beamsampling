# Calculate residuals for models

import data_handling as dh
import pickle
from os.path import join as pj
import os
from definitions import OUTPUT_DIR
import argparse
import utils
import ideal_observer as i_o
import linear_markov
import numpy as np

models = {'iHMM': i_o.IdealObserverSamples,
          'LinearMarkov': linear_markov.LinearMarkovSamples,
          'Markov': i_o.IdealObserverSamples}

def residuals(participant, ini, epoch_train, epoch_test, commit, model_class, output = None):
    data = dh.import_data(participant=participant, epoch=epoch_test)
    model = model_class(participant=participant, ini=ini, commit=commit, epoch=epoch_train)
    data['rt'] = model.predict_rt(data['Y'])
    
    if output is not None:
        with open(output,'wb') as file:
            pickle.dump(data, file)

    return(data)


def residuals_into_folder(participants, inis, epoch_trains, epoch_tests, commit, model_class, keep_mean=False):
    folder = pj(OUTPUT_DIR,'residuals_'+commit+'_'+model_class.__name__)
    if not os.path.exists(folder):
        os.makedirs(folder)
        
    for epoch_train in epoch_trains:
        prefix = 'residual_'+model_class.__name__+'_'+epoch_train+'_'+commit+'_'+('km_' if keep_mean else '')
        print('Epoch (train) ', epoch_train,':   [', end='')
        for participant in participants:
            for epoch_test in epoch_tests:
                input_filename = dh.data_filename(participant, epoch_test)
                data = dh.import_data(participant=participant, block=epoch_test)
                for ini in inis:
                    kwargs = dict(last_n_samples=60)
                    if model_class == 'Markov':
                        ini = ini+'_markov'
                    model = model_class(participant=participant, ini=ini, commit=commit, epoch=epoch_train, **kwargs)
                    out = data.copy()
                    out['rt'] = data['rt']-model.predict_rt(data['Y'])
                    if keep_mean:
                        out['rt'] = out['rt'] + np.mean(data['rt']) - np.mean(out['rt'])
                    with open(pj(folder, prefix+input_filename),'wb') as file:
                        pickle.dump(out, file)
                print('.', end='')
        print(']')
        
def main():
    parser = argparse.ArgumentParser(description='Compute residuals given model samples.')
    parser.add_argument('-e', metavar = 'experiment', type=str,
                        help='Experiment')
    parser.add_argument('-p', metavar = 'participants', type=str,
                        help='Participants (comma separated)')
    parser.add_argument('-i', metavar='ini_files', type=str,
                        help='Ini files (of the models used), comma separated.')
    parser.add_argument('--e_trains', metavar='epoch_trains', type=str,
                        help='Epochs in which training is taken place (comma separated).')
    parser.add_argument('--e_tests', metavar = 'epoch_tests', type=str,
                        help='Epochs in which calculate residuals (comma separated).')
    parser.add_argument('-c', metavar = 'commit', type=str,
                        help='Commit (latest allowed)')
    parser.add_argument('-m', metavar = 'model_class', type=str,
                        help='Model class: iHMM, LinearMarkov, Markov')
    parser.add_argument('--keep_mean', action='store_true',
                        help='Keep mean of RT intact.')
    
    args = parser.parse_args()
    participants = args.p.split(',')
    experiment = args.e
    participants = [experiment+'_'+p for p in participants]
    inis = args.i.split(',')
    inis = [ini+'.ini' for ini in inis]
    epoch_trains = args.e_trains.split(',')
    epoch_tests = args.e_tests.split(',')
    commit = args.c
    if commit == 'latest':
        commit = utils.get_git_hash('..')[:7]
    model_class = args.m
    print(participants)
    print(inis)
    print(epoch_trains)
    residuals_into_folder(participants, inis, epoch_trains, epoch_tests, commit, models[model_class], args.keep_mean)
        
if __name__=='__main__':
    main()