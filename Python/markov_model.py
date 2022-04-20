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

sm_pi_rtparams = pickle.load(open(utils.find_file('sm_pi_rtparams.pkl'),'rb'))
sm_linear_markov = pickle.load(open(utils.find_file('sm_linear_markov.pkl'),'rb'))

LATER_MODEL = 'LATER_MODEL'
LINEAR_MODEL = 'LINEAR_MODEL'

def sample_markov(Y, 
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
                     K = np.max(Y)+1,
                     D = np.max(Y)+1,
                     T = len(Y),
                     tau0_scale = hypers['tau0_scale'],
                     tau0_shape = hypers['tau0_shape'],
                     mu_scale = hypers['mu_scale'],
                     mu_shape = hypers['mu_shape'],
                     sigma_scale = hypers['sigma_scale'],
                     sigma_shape = hypers['sigma_shape'],
                     error_sigma = hypers['error_sigma'],
                     rt_max = hypers['rt_max'],
                     alpha0 = hypers['alpha0'],
                     phi = np.eye(np.max(Y)+1),
                     rt = rt,
                     correct_response = correct_response,
                     filters = np.int8(_filters),
    )
    time_stan_start = time.time()
    fit = sm_pi_rtparams.sampling(data = stan_data,
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
    
    for tau0, mu, sigma, pi, jll in zip(S['tau0'], S['mu'], S['sigma'], S['pi'], S['lp__']):
        sample = dict()
        sample['tau0'] = tau0
        sample['mu'] = mu
        sample['sigma'] = sigma
        sample['pi'] = pi
        sample['phi'] = np.eye(np.max(Y)+1)
        sample['jll'] = jll
        samples.append(sample)
        stats['jll'].append(jll)
        stats['stan_time'].append(time_stan_finish-time_stan_start)
    if output is not None:
        with open(output,'wb') as file:
                pickle.dump(dict(stats = stats,
                                 samples = samples, 
                                 hypers = hypers,
                                 sampling_params = sampling_params,
                                 model_type=LATER_MODEL),
                            file)
                file.close()
    return(samples)

def sample_linear_markov(Y, rt, correct_response, _filters, hypers, sampling_params, output = None):
    stan_data = dict(Y = Y+1,
                     D = np.max(Y)+1,
                     T = len(Y),
                     mu_mean = hypers['mu_mean'],
                     mu_sd = hypers['mu_sd'],
                     sigma_mean = hypers['sigma_mean'],
                     sigma_sd = hypers['sigma_sd'],
                     pi_sd = hypers['pi_sd'],
                     rt = rt,
                     correct_response = correct_response,
                     filters = np.int8(_filters),
    )
    time_stan_start = time.time()
    fit = sm_linear_markov.sampling(data = stan_data,
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
    
    for mu, sigma, pi_matrix, jll in zip(S['mu'], S['sigma'], S['pi_matrix'], S['lp__']):
        sample = dict()
        sample['mu'] = mu
        sample['sigma'] = sigma
        sample['pi'] = pi_matrix
        sample['phi'] = np.eye(np.max(Y)+1)
        sample['jll'] = jll
        samples.append(sample)
        stats['jll'].append(jll)
        stats['stan_time'].append(time_stan_finish-time_stan_start)
    if output is not None:
        with open(output,'wb') as file:
                pickle.dump(dict(stats = stats,
                                 samples = samples, 
                                 hypers = hypers,
                                 sampling_params = sampling_params,
                                 model_type=LINEAR_MODEL),
                            file)
                file.close()
    return(samples)

def sample_markov_from_file(config_filename, data_filename = None, model_type=LATER_MODEL):
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
    print('------       MARKOV MODEL INFERENCE       -------')
    print('-------------------------------------------------')
    print('MODEL TYPE: ',model_type)
    print('Initiating sampling using config file:', config_filename)
    
    if data_filename is not None:
        with open(os.path.join(DATA_DIR,data_filename),'rb') as file:
            data = pickle.load(file)
        
        if model_type == LATER_MODEL:
            output_file = os.path.join(OUTPUT_DIR,utils.get_git_hash('..')[:7]+'_'+config_filename[:-4]+'_'+data_filename[:-4]+'_markov_samples.pkl')
        if model_type == LINEAR_MODEL:
            output_file = os.path.join(OUTPUT_DIR,utils.get_git_hash('..')[:7]+'_'+config_filename[:-4]+'_'+data_filename[:-4]+'_linear_markov_samples.pkl')
        print('Importing stimulus and reaction time data from', data_filename)
        print('Writing results into', output_file)
    
    else:
        print('No data given.')
    
    print('# Hyperparameter settings #')
    hypers = dict()
    for key in config['HYPERS']:
        hypers[key] = np.float32(config['HYPERS'][key])    
        print('    ',key,'=',hypers[key])
    print()
    
    sampling_params = dict()
    print('# Sampling parameters #')
    for key in config['SAMPLING_PARAMS']:
        sampling_params[key] = int(config['SAMPLING_PARAMS'][key])    
        print('    ',key,'=',sampling_params[key])
    print()
    
    rt = data['rt']
    Y = data['Y']
    
    if 'correct_response' in data:
        correct_response = data['correct_response']
    else:
        correct_response = np.ones(len(rt), dtype = np.int32)

    _filters = filters.get_filter(data,config)
    
    if model_type == LATER_MODEL:
        sample_markov(Y, 
                      rt, 
                      correct_response, 
                      _filters, 
                      hypers, 
                      sampling_params, 
                      output = output_file)
        
    if model_type == LINEAR_MODEL:
        sample_linear_markov(Y, rt, correct_response, _filters, hypers, sampling_params, output = output_file)
        
    print('--------------------------')
    print('--- SAMPLING FINISHED. ---')
    print('--------------------------')
    print()