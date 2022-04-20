## Wrapper to infer internal models that are purely Markov
## a.k.a. models where Phi is unity and K = len(unique(Y))

import pickle
import pystan
import numpy as np
import utils
import configparser
import os
import time
from definitions import DATA_DIR, CONFIG_DIR, OUTPUT_DIR
import filters

sm_hmm = pickle.load(open(utils.find_file('sm_hmm.pkl'),'rb'))

def sample_hmm(Y,
               rt,
               correct_response,
               _filters,
               hypers,
               sampling_params,
               output = None):

    for i in range(len(rt)):
        if rt[i] < 1:
            rt[i] = 1

    correct_response = correct_response * (np.array(rt)>1)

    stan_data = dict(Y = Y+1,
                     K = hypers['K'],
                     D = np.max(Y)+1,
                     T = len(Y),
                     M = np.zeros((hypers['K'], np.max(Y)+1)),
                     H = hypers['H'],
                     tau0_scale = hypers['tau0_scale'],
                     tau0_shape = hypers['tau0_shape'],
                     mu_scale = hypers['mu_scale'],
                     mu_shape = hypers['mu_shape'],
                     sigma_scale = hypers['sigma_scale'],
                     sigma_shape = hypers['sigma_shape'],
                     error_sigma = hypers['error_sigma'],
                     rt_max = hypers['rt_max'],
                     alpha0 = hypers['alpha0'],
                     rt = rt,
                     correct_response = correct_response,
                     filters = np.int8(_filters),
                     new_stream = hypers['new_stream'],
    )
    time_stan_start = time.time()
    fit = sm_hmm.sampling(data = stan_data,
                          iter = sampling_params['nums']+sampling_params['numb'],
                          seed = [int(sampling_params['seed'])],
                          #warmup = sampling_params['numb'],
                          chains = 1,
                         )
    time_stan_finish = time.time()

    S = fit.extract()
    samples = []
    stats = dict(jll = [],
                 stan_time=[])

    for tau0, mu, sigma, pi, phi, jll in zip(S['tau0'], S['mu'], S['sigma'], S['pi'], S['phi'], S['lp__']):
        sample = dict()
        sample['tau0'] = tau0
        sample['mu'] = mu
        sample['sigma'] = sigma
        sample['pi'] = pi
        sample['phi'] = phi
        sample['jll'] = jll
        samples.append(sample)
        stats['jll'].append(jll)
        stats['stan_time'].append(time_stan_finish-time_stan_start)
    if output is not None:
        with open(output,'wb') as file:
            pickle.dump(dict(stats = stats,
                             samples = samples,
                             hypers = hypers,
                             sampling_params = sampling_params),
                        file)
            file.close()
    return(samples)

def sample_hmm_from_file(config_filename, data_filename, K):
    print(config_filename)
    config = configparser.ConfigParser()
    config.read(os.path.join(CONFIG_DIR,config_filename))

    print('-------------------------------------------------')
    print('--- iHMM Cognitive Tomography by Balazs Torok ---')
    print('-------------------------------------------------')
    print()
    print('% Thanks to J Van Gael for providing code of iHMM Beam Sampling')
    print('% This code makes heavy use of his work')
    print()
    print('-------------------------------------------------')
    print('--------      HMM MODEL INFERENCE        --------')
    print('-------------------------------------------------')
    print('NUMBER OF STATES: ',K)
    print('Initiating sampling using config file:', config_filename)

    if data_filename is not None:
        with open(os.path.join(DATA_DIR,data_filename),'rb') as file:
            data = pickle.load(file)

        output_file = os.path.join(OUTPUT_DIR,utils.get_git_hash('..')[:7]+'_'+config_filename[:-4]+'_'+data_filename[:-4]+'_hmm'+str(K)+'_samples.pkl')
        print('Importing stimulus and reaction time data from', data_filename)
        print('Writing results into', output_file)

    else:
        print('No data given.')

    print('# Hyperparameter settings #')
    hypers = dict(K=K)
    for key in config['HYPERS']:
        hypers[key] = np.float32(config['HYPERS'][key])
        print('    ',key,'=',hypers[key])
    print()

    if 'new_stream_at_lowest_trial_in_block' not in hypers:
        hypers['new_stream_at_lowest_trial_in_block'] = 0

    if hypers['new_stream_at_lowest_trial_in_block'] == 1:
        hypers['new_stream'] = filters.new_stream(data)
    else:
        hypers['new_stream'] = np.zeros(len(data['rt']),dtype=np.int8)
        hypers['new_stream'][0] = 1

    print(hypers['new_stream'])

    sampling_params = dict()
    print('# Sampling parameters #')
    for key in config['SAMPLING_PARAMS']:
        sampling_params[key] = int(config['SAMPLING_PARAMS'][key])
        print('    ',key,'=',sampling_params[key])
    print()

    rt = data['rt']
    Y = data['Y']
    hypers['H'] = np.array([np.repeat(1,np.max(Y)+1)])*hypers['h']
    if 'correct_response' in data:
        correct_response = data['correct_response']
    else:
        correct_response = np.ones(len(rt), dtype = np.int32)

    _filters = filters.get_filter(data,config)

    sample_hmm(Y,
               rt,
               correct_response,
               _filters,
               hypers,
               sampling_params,
               output = output_file)

    print('--------------------------')
    print('--- SAMPLING FINISHED. ---')
    print('--------------------------')
    print()
