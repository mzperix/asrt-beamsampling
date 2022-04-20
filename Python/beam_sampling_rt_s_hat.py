### BEAM SAMPLING FUNCTIONS

### CONDITIONING ON RT !!!!


### Based on: Matlab beamsampling functions

import numpy as np
import scipy.stats as sp
import sys, os
from scipy.special import gammaln
import carpenter_williams_model as cwm
import pickle
import pystan
import time
import utils
from wurlitzer import sys_pipes
import configparser
from ideal_observer import IdealObserver
import filters
from definitions import DATA_DIR, OUTPUT_DIR, CONFIG_DIR

# TEMPORARY FOR DEBUGGING
import matplotlib.pyplot as plt

#RANDOM_SAMPLER = 'MATLAB'
RANDOM_SAMPLER = 'NUMPY'

if RANDOM_SAMPLER == 'MATLAB':
    print('Using Matlab samplers.')
    import matlab_sampler as ms
    random_gamma = ms.gamrnd
    random_beta  = ms.betarnd
    random_uniform = ms.rand
    random_binomial = ms.binornd

if RANDOM_SAMPLER == 'NUMPY':
    random_gamma = np.random.gamma
    random_beta = np.random.beta
    random_uniform = np.random.uniform
    random_binomial = np.random.binomial

# LOAD STAN MODELS
#sm_emission_matrix = pickle.load(open(utils.find_file('sm_emission_matrix.pkl'),'rb'))
#sm_rt_model = pickle.load(open(utils.find_file('sm_rt_model.pkl'),'rb'))
#sm_joint_e_rt = pickle.load(open(utils.find_file('sm_joint_e_rt.pkl'),'rb'))
sm_joint_pi_phi_rt = pickle.load(open(utils.find_file('sm_pi_phi_rtparams.pkl'),'rb'))
sm_hypers = pickle.load(open(utils.find_file('sm_hypers.pkl'),'rb'))

def check_error(sample):
    if sample['tau0'] > 1e2:
        return True
    if sample['tau0'] < 1e-4:
        return True
    if sample['mu'] > 1e2:
        return True
    if sample['mu'] < 1e-5:
        return True
    if sample['sigma'] > 1e2:
        return True
    if sample['sigma'] < 1e-9:
        return True
    else:
        return False

def sample_rt_parameters(S, Y, hypers, rt, rt_params, phi, numi_rt):
    T = np.size(S)
    log_pred_prob = np.zeros(T)

    log_phi = np.log(1e-15+phi)

    for t in range(T):
        log_pred_prob[t] = log_phi[S[t],Y[t]]

    data = dict(N = T,
                rt = rt,
                logp = log_pred_prob,
                tau0_shape = hypers['tau0_a'],
                tau0_scale = 1.0 / hypers['tau0_b'],
                mu_shape = hypers['mu_a'],
                mu_scale = 1.0 / hypers['mu_b'],
                sigma_shape = hypers['sigma_a'],
                sigma_scale = hypers['sigma_b'],
                error_sigma = 1,
                rt_max = hypers['rt_max'])

    fit = sm_rt_model.sampling(data = data, 
                           init = [dict(tau0 = rt_params['tau0'],
                                       mu = rt_params['mu'],
                                       sigma = rt_params['sigma'])],
                           chains = 1, iter = numi_rt, warmup = 1, 
                           verbose = False)
    samples = fit.extract()
    sample_rt_params = dict(tau0 = samples['tau0'][-1],
                            mu = samples['mu'][-1],
                            sigma = samples['sigma'][-1])
    return(sample_rt_params)


# Matlab function: dirichlet_sample.m
def sample_dirichlet(alpha):
    sample = random_gamma(alpha, 1.0)
    sample = sample / np.sum(sample)
    return(sample)


# Matlab function: SampleEmissionMatrix(S, Y, K, H)
def sample_emission_matrix(S, Y, K, H, rt, rt_params, phi0, numi_phi):
    # S: state sequence (np.array)
    # Y: emission sequence (np.array)
    # K: number of represented states
    # H: parameter of emission prior (~ Dir(H))
    # rt: reaction times (np.array)
    # rt_params: parameters of rt model (CW)

    L = np.size(H, 1)
    M = np.zeros((K,L))
    
    for k in range(K):
        for l in range(L):
            M[k,l] = np.sum(Y[S==k]==l)
    
    data = dict(T = np.size(S),
                D = L,
                K = K,
                S = S+1, # array representation in stan goes from 1, and we look at Phi[S[t],Y[t]]
                Y = Y+1,
                rt = rt,
                tau0 = rt_params['tau0'],
                mu = rt_params['mu'],
                sigma = rt_params['sigma'],
                M = M,
                H = H)

    fit = sm_emission_matrix.sampling(data=data, iter=numi_phi, chains=1, warmup = 1, verbose=True)
    Phi = fit.extract()['phi'][-1]
    return(Phi)


def sample_joint_e_rt(S, Y, K, H, rt, hypers, numi, phi0, rt_params):
    # Sample emission and rt parameters TOGETHER (using stan)
    L = np.size(H, 1)
    M = np.zeros((K,L))
    
    for k in range(K):
        for l in range(L):
            M[k,l] = np.sum(Y[S==k]==l)

    data = dict(N = np.size(S),
                rt = rt,
                
                tau0_shape = hypers['tau0_a'],
                tau0_scale = 1.0 / hypers['tau0_b'],
                mu_shape = hypers['mu_a'],
                mu_scale = 1.0 / hypers['mu_b'],
                sigma_shape = hypers['sigma_a'],
                sigma_scale = hypers['sigma_b'],
                error_sigma = 1,
                rt_max = hypers['rt_max'],

                T = np.size(S),
                D = L,
                K = K,
                S = S+1, # array representation in stan goes from 1, and we look at Phi[S[t],Y[t]]
                Y = Y+1,
                
                M = M,
                H = H)

    fit = sm_joint_e_rt.sampling(data=data, 
                                 iter = numi, 
                                 chains = 1, 
                                 warmup = int(numi / 2), 
                                 init = [dict(phi=phi0, 
                                              tau0 = rt_params['tau0'],
                                              mu = rt_params['mu'],
                                              sigma = rt_params['sigma'])],
                                 verbose = True)
    samples = fit.extract()
    Phi = samples['phi'][-1]
    tau0 = samples['tau0'][-1]#np.mean(samples['tau0'])
    mu = samples['mu'][-1]#np.mean(samples['mu'])
    sigma = samples['sigma'][-1]#np.mean(samples['sigma'])
    return(tau0, mu, sigma, Phi)


def sample_pi_phi_rtparams(Y, K, H, rt, correct_response, _filters, epsilon, sample, hypers, numi, init):
    # Sample pi, phi, rt parameters given slicing variables and hyper parameters
    # Y: discrete stimuli
    # H: prior for emission matrix
    # epsilon: slicing parameters
    # hypers: hyper parameters
    #   alpha: 
    #   beta:
    #   gamma:
    #   tau0_scale, tau0_shape:   prior parameters for tau0
    #   mu_shape, mu_scale:       prior parameters for mu
    #   sigma_shape, sigma_scale: prior parameters for sigma
    #   error_sigma: measurement error sd
    #   beta_softmax: parameter of softplus and softmax functions for pi
    # rt_max: maximal recorded RT (cutoff of RT distribution)
    # rt: reaction times
    # K: number of states currently taken into consideration = (dimension of beta)-1
     
    D = np.max(Y)+1
    T = np.size(Y)
    
    stan_data = dict(T = T,
                     K = K,
                     tau0_scale = hypers['tau0_scale'],
                     tau0_shape = hypers['tau0_shape'],
                     mu_shape = hypers['mu_shape'],
                     mu_scale = hypers['mu_scale'],
                     sigma_shape = hypers['sigma_shape'],
                     sigma_scale = hypers['sigma_scale'],
                     error_sigma = hypers['error_sigma'],
                     rt_max = hypers['rt_max'],
                     rt = rt,
                     correct_response = correct_response,
                     filters = _filters,
                     new_stream = hypers['new_stream'],
                     D = D,
                     Y = Y+1,
                     epsilon = epsilon,
                     H = hypers['H'],
                     M = np.zeros((K, np.max(Y)+1)),
                     beta_softmax = hypers['beta_softmax'],
                     alpha0 = sample['alpha0'],
                     beta   = np.squeeze(sample['Beta']),
                     )
    
    fit = sm_joint_pi_phi_rt.sampling(data = stan_data, 
                                      iter = numi, 
                                      chains = 1, 
                                      #warmup = 1, 
                                      init = [dict(pi = init['pi'],
                                                   phi = init['phi'], 
                                                   tau0 = init['tau0'],
                                                   mu = init['mu'],
                                                   sigma = init['sigma'])],
                                      verbose = True,
                                      seed = np.random.get_state()[2])
    samples = fit.extract()
    Pi = samples['pi'][-1]
    Phi = samples['phi'][-1]
    tau0 = samples['tau0'][-1]#np.mean(samples['tau0'])
    mu = samples['mu'][-1]#np.mean(samples['mu'])
    sigma = samples['sigma'][-1]#np.mean(samples['sigma'])
    return(Pi, Phi, tau0, mu, sigma)


# Matlab function: SampleTransitionMatrix(S, H)
def sample_transition_matrix(S, H):
    # S: state sequence (np.array)
    # H: parameter of dirichlet prior

    K = np.size(H,1)
    T = int(np.size(S))
    Pi = np.zeros((K,K))
    N = np.zeros((K,K))
    for t in range(1,T):
        N[S[t-1], S[t]] += 1

    for k in range(K):
        Pi[k,:] = sample_dirichlet(N[k,:]+H)
    return(Pi)


# Matlab function: iHmmJointLogLikelihood (S, Y, Beta, alpha0, H)
def ihmm_joint_llike(Y, rt, correct_response, _filters, Pi, Phi, tau0, mu, sigma, alpha0, gamma, Beta, epsilon, hypers):
    # Y: emission sequence (np.array)
    # rt: reaction times (np.array)
    # Pi: transition matrix (np.array)
    # Pi: emission matrix (np.array)
    # tau0: rt parameter
    # mu: rt parameter
    # sigma: rt parameter
    # alpha0: hyperparameter
    # gamma: hyperparameter
    # Beta: parameter of transition matrix prior (result of DP)
    # hypers: hyperparameters defining the priors
    
    K = np.size(Pi,0)
    T = np.size(Y)

    # Compute the log likelihood
    logp = 0
    
    # Hyperparameters (SKIPPED ATM.)
    #logp += sp.gamma.logpdf(tau0, hypers['alpha0_a'], scale = hypers['alpha0_scale'])
    #logp += sp.dirichlet.logpdf(Beta, alpha)
    
    # Prior for reaction time parameters
    logp += sp.gamma.logpdf(tau0, hypers['tau0_shape'], scale = hypers['tau0_scale'])
    logp += sp.gamma.logpdf(tau0, hypers['mu_shape'], scale = hypers['mu_scale'])
    logp += sp.gamma.logpdf(tau0, hypers['sigma_shape'], scale = hypers['sigma_scale'])
    
    # Reaction times
    observer = IdealObserver(tau0 = tau0,
                             mu = mu,
                             sigma = sigma,
                             Pi = Pi,
                             Phi = Phi)
    logp += np.sum(correct_response*_filters*observer.loglike(Y, rt))

    # Pi (transition matrix)
    for k in range(K):
        #logp += sp.dirichlet.logpdf(Pi[k,:], alpha0*Beta[0,:])
        logp += sp.dirichlet.logpdf(Pi[k,:], alpha0*np.concatenate([np.repeat(1/K,K),np.array([epsilon])]))
        
    # Phi (emission matrix)
    for k in range(K):
        Phi[k,:] += 1e-3
        Phi[k,:] /= np.sum(Phi[k,:])
        try:
            logp += sp.dirichlet.logpdf(Phi[k,:], hypers['H'][0,:])
        except:
            print('Phi: ', Phi[k,:])
            print('H: ', hypers['H'][0,:])

    return(logp)


# Matlab function: iHmmHyperSample
# Resamples hyperparameters of an infinite hmm
def sample_ihmm_hyper(hypers, numi, pi = None, K = None):
#   [beta, alpha0, gamma] = ...
#   iHmmHyperSample(S, ibeta, ialpha0, igamma, hypers, numi) resamples the
#   hyperparameters given the state sequence S, the previous
#   hyperparameters ibeta, ialpha0, igamma and their respective
#   hyper-hyperparameters in the structure hypers (needs alpha0_a,
#   alpha0_b, gamma_a and gamma_b fields corresponding to gamma prior on
#   the hyperparameters). If the hyper-hyperparameters are not given,  the
#   estimated alpha0 and gamma will be the same as the input alpha0 and
#   gamma. numi is the number of times we run the Gibbs samplers for alpha0
#   and gamma (see HDP paper or Escobar & West); we recommend a value of
#   around 20. The function returns the new hyperparameters, the CRF counts
#   (N) and the sampled number of tables in every restaurant (M).
#
#   Note that the size of the resampled beta will be the same as the size
#   of the original beta.
    if pi is None:
        # Sample from prior
        if K is None:
            raise RuntimeError('Nor pi, nor K is given')
        else:
            if 'alpha0' in hypers:
                alpha0 = hypers['alpha0']
            else:
                alpha0 = random_gamma(hypers['alpha0_a'], hypers['alpha0_b'])
            if 'gamma' in hypers:
                gamma = hypers['gamma']
            else:
                gamma  = random_gamma(hypers['gamma_a'], hypers['gamma_b'])
            beta = np.ones((1,K+1))
            for k in range(1,(K+1)):
                be = beta[0,k-1]
                bg = random_beta(1.0, gamma)
                beta[0,k-1] = bg * be
                beta[0,k] = (1.0-bg) * be
           
    else:
        K = np.size(pi,0)  # number of states in iHmm
        fit = sm_hypers.sampling(data = dict(K = K,
                                             alpha0_a = hypers['alpha0_a'],
                                             alpha0_b = hypers['alpha0_b'],
                                             gamma_a = hypers['gamma_a'],
                                             gamma_b = hypers['gamma_b'],
                                             pi = pi),
                                iter = numi,
                                chains = 1,
                                seed = np.random.get_state()[2])
        samples = fit.extract()
        if 'alpha0' in hypers:
            alpha0 = hypers['alpha0']
        else:
            alpha0 = samples['alpha0'][-1]
        if 'gamma' in hypers:
            gamma = hypers['gamma']
        else:
            gamma  = samples['gamma'][-1]
        beta = np.reshape(samples['beta'][-1], (1,K+1))
    
    return(beta, alpha0, gamma)


def sample_beam_ihmm(K0, Y, rt, correct_response, _filters, hypers, sampling_params, output = None, with_plots = False):
# % sample_beam_ihmm samples internal models and rt parameters from 
# % rt measurements and inputs (Y)
# %
# %   Input Parameters:
# %   - K0: initial number of states
# %   - Y: training sequence of arbitrary length,
# %   - rt: array of reaction times
# %   - correct_response: array of correct responses
# %   - _filters: array of other filters (0: leave out from conditioning, 1: use in conditioning)
# %   - hypers: a structure that describes the hyperparameters for the sampler.
# %             REQUIRED:
# %                 - alpha0 prior parameters: alpha0_a, alpha0_b (gamma dist.)
# %                 - gamma  prior parameters: gamma_a, gamma_b (gamma dist.)
# %                 - rt parameter priors: tau0_shape, tau0_scale: gamma dist.
# %                                        mu_shape, mu_scale: gamma dist.
# %                                        sigma_shape, sigma_scale: gamma dist.
# %             OPTIONAL:
# %                 - epsilon0: the initial slicing parameter value
# %   - sampling_params:
# %             - numb: number of burnin iterations
# %             - numi_pi_phi_rtparams: number of STAN steps for phi_pi_rtparams sampling
# %             - numi_hypers: number of STAN steps for hyper sampling
# %             - nums: number of samples to output
# %
# %   Output Parameters:
# %   - S: is a list of sample dicts, including: Pi, Phi, tau0, mu, sigma, K
# %   - stats: is a structure that contains a variety of statistics for every
# %            iteration of the sampler: K, alpha0, gamma, the size of the
# %            trellis and the marginal likelihood.

    np.random.seed(sampling_params['seed'])

    if output is None:
        output = input('Output filename:')

    # Initialize the sampler.
    T = np.size(Y) # Number of timesteps

    # Setup structures to store the output
    S = []
    stats = dict(git_commit_hash = utils.get_git_hash('..'),
                 K = np.zeros(sampling_params['numb'] + sampling_params['nums']),
                 alpha0 = np.zeros(sampling_params['numb'] + sampling_params['nums']),
                 gamma  = np.zeros(sampling_params['numb'] + sampling_params['nums']),
                 jll = np.zeros(sampling_params['numb'] + sampling_params['nums']),
                 stan_time = np.zeros(sampling_params['numb'] + sampling_params['nums']))

    # Initialize hypers
    sample = dict()
    sample['K'] = K0
    sample['Beta'], sample['alpha0'], sample['gamma'] = sample_ihmm_hyper(hypers = hypers,
                                                                          numi = sampling_params['numi_hypers'],
                                                                          K = K0, 
                                                                         )
    
    # Initialize rt parameters
    sample['tau0'] = random_gamma(hypers['tau0_shape']) * hypers['tau0_scale']
    sample['mu'] = np.mean((sample['tau0']-np.log(0.25)) / rt)
    sample['sigma'] = 40*np.std(((sample['tau0']-np.log(0.25)) / rt))
    print('RT params initialized as:  tau0:', sample['tau0'], ' mu:', sample['mu'], ' sigma:', sample['sigma'])
    print('Hyperparams initialized as: alpha0:', sample['alpha0'], ' gamma:', sample['gamma'])
    
    # Sample the emission and transition probabilities.
    sample['Pi'] = np.vstack([np.array(sample_dirichlet(sample['Beta']*sample['alpha0'])) for i in range(sample['K'])])
#   sample['Pi'] = np.hstack([np.ones((sample['K'],sample['K']))/sample['K'], np.zeros((sample['K'],1))])
    sample['Phi'] = np.vstack([np.array(sample_dirichlet(hypers['H'])) for i in range(sample['K'])])
#   sample['Phi'] = np.ones((sample['K'], 4)) / 4
    
    if with_plots:
        print('Pi and Phi')
        fig, axes = plt.subplots(1,2)
        axes[0].matshow(sample['Pi'], vmin = 0, vmax = 1)
        axes[1].matshow(sample['Phi'], vmin = 0, vmax = 1)
        plt.show()

    # Set up epsilon
    if 'epsilon0' in hypers:
        epsilon = hypers['epsilon0']
    else:
        epsilon = 0.1
    
    iteration = 0
    while iteration <= (sampling_params['numb'] + sampling_params['nums'] - 1):
        # Safety check
        assert(np.size(sample['Phi'],0) == np.size(sample['Beta'],1) - 1)

        # Sample auxilary variable epsilon
        epsilon = (random_uniform()
                   *(hypers['epsilon_max']-hypers['epsilon_min'])
                   +hypers['epsilon_min']) #2e-2+(random_uniform() * np.unique(sample['Pi'] * np.array(sample['Pi'] > epsilon))[1])
        
        while np.max(sample['Pi'][:, -1]) > epsilon: # Break the Pi[k] stick some more.
            pl = np.size(sample['Pi'],1)
            bl = np.size(sample['Beta'],1)

            # Safety check
            assert(bl == pl)

            # Add row to transition matrix.
            sample['Pi'] = np.vstack([sample['Pi'], sample_dirichlet(sample['alpha0'] * sample['Beta'])])
            sample['Phi']= np.vstack([sample['Phi'], sample_dirichlet(hypers['H'])])

            # Break the beta stick
            be = sample['Beta'][0,-1]
            bg = random_beta(1.0, sample['gamma'])
            sample['Beta'][0,bl-1] = bg * be;
            sample['Beta'] = np.hstack([sample['Beta'], np.array([[(1.0-bg) * be]])])
            
            pe = sample['Pi'][:,-1].copy()
            a = np.tile(sample['alpha0']*sample['Beta'][0,-2], (bl, 1))
            b = sample['alpha0'] * (1.0 - np.sum(sample['Beta'])+sample['Beta'][0,-1])
            if (np.min(a) < 1e-2) or (np.min(b) < 1e-2):
                pg = np.transpose(random_binomial(1.0, a / (a+b)))[0]
            else:
                pg = random_beta(a, b)
            pg = np.array([np.float64(x) for x in pg])
            sample['Pi'][:, pl-1] = pg * pe
            sample['Pi'] = np.hstack([sample['Pi'], np.transpose(np.array([(1.0-pg) * pe], dtype = np.float32))]
                                     )

        sample['K'] = np.size(sample['Pi'],0)
        
        sample['Pi'] += 1e-4
        
        ## TEST PI
        for k in range(sample['K']):
            sample['Pi'][k,:] = np.divide(sample['Pi'][k,:],np.sum(sample['Pi'][k,:]))
            sample['Pi'][k,-1] = 1-np.sum(sample['Pi'][k,:-1])
            #print('%.30f' % (np.sum(sample['Pi'][k,:])-1))
            #print(sample['Pi'][k,:])
            assert np.allclose(np.sum(sample['Pi'][k,:]), 1.0, atol = 1e-7, rtol = 1e-7)
        
        if with_plots:
            print('Pi and Phi after stick breaking')
            fig, axes = plt.subplots(1,2)
            axes[0].matshow(sample['Pi'], vmin = 0, vmax = 1)
            axes[1].matshow(sample['Phi'], vmin = 0, vmax = 1)
            plt.show()
            print('And Beta')
            plt.plot(np.squeeze(sample['Beta']), '.')
            plt.show()

        # Cleanup our state space by removing redundant states.
        cleanup = True
        while cleanup:
            i = np.argmin(np.max(sample['Pi'][:,:-1], axis = 0))
            # If we can put the lowest goto probability column (i.e. state)
            # into the 'unrepresented' states (last column), we delete it
            if np.max(sample['Pi'][:,-1]+sample['Pi'][:,i]) < epsilon:
                # Delete state i
                sample['Beta'][0,-1] += sample['Beta'][0,i]
                sample['Beta'] = np.delete(sample['Beta'], i, axis = 1)
                sample['Pi'][:,-1] += sample['Pi'][:,i]
                sample['Pi'] = np.delete(sample['Pi'], i, axis = 0)
                sample['Pi'] = np.delete(sample['Pi'], i, axis = 1)
                sample['Phi'] = np.delete(sample['Phi'], i, axis = 0)
                sample['K'] = np.size(sample['Pi'],0)
            else:
                cleanup = False
                
        if with_plots:
            print('Pi and Phi after cleaning up')
            fig, axes = plt.subplots(1,2)
            axes[0].matshow(sample['Pi'], vmin = 0, vmax = 1)
            axes[1].matshow(sample['Phi'], vmin = 0, vmax = 1)
            plt.show()
            print('And Beta')
            plt.plot(np.squeeze(sample['Beta']), '.')
            plt.show()
        
        # Safety check
        assert(sample['K'] == np.size(sample['Beta'])-1)
        assert(sample['K'] == np.size(sample['Phi'],0))

        # Resample Pi, Phi and RT params
        # print('loglike before sampling Phi and rt params: ', ihmm_joint_llike(sample['S'], Y, sample['Beta'], sample['alpha0'], hypers['H'], rt_params, rt, sample['Phi']))
        time_stan_start = time.time()
        sample['Pi'], sample['Phi'], sample['tau0'], sample['mu'], sample['sigma'] = sample_pi_phi_rtparams(Y = Y, 
                                                                                                            K = sample['K'], 
                                                                                                            H = hypers['H'], 
                                                                                                            rt = rt, 
                                                                                                            correct_response = correct_response,
                                                                                                            _filters = np.int8(_filters), 
                                                                                                            epsilon = epsilon,
                                                                                                            sample = sample,
                                                                                                            hypers = hypers, 
                                                                                                            init = dict(pi = sample['Pi'], 
                                                                                                                        phi = sample['Phi'],
                                                                                                                        tau0 = sample['tau0'],
                                                                                                                        mu = sample['mu'],
                                                                                                                        sigma = sample['sigma'],
                                                                                                                       ),
                                                                                                            numi = sampling_params['numi_pi_phi_rtparams'])
        # print('loglike after  sampling Phi and rt params: ', ihmm_joint_llike(sample['S'], Y, sample['Beta'], sample['alpha0'], hypers['H'], rt_params, rt, sample['Phi']))
        time_stan_finish = time.time()
        if with_plots:
            print('Pi and Phi after STAN')
            fig, axes = plt.subplots(1,2)
            axes[0].matshow(sample['Pi'], vmin = 0, vmax = 1)
            axes[1].matshow(sample['Phi'], vmin = 0, vmax = 1)
            plt.show()
        
        # Resample Beta given the transition probabilities.
        sample['Beta'], sample['alpha0'], sample['gamma'] =  sample_ihmm_hyper(hypers = hypers,
                                                                               numi = sampling_params['numi_hypers'],
                                                                               pi = sample['Pi'], 
                                                                              )
        # Safety checks.
        assert(np.size(sample['Pi'],0) == sample['K'])
        assert(np.size(sample['Pi'],1) == sample['K']+1)
        assert(sample['K'] == np.size(sample['Beta'])-1)
        assert(np.min(sample['Pi'] >= 0))

        # Prepare next iteration.
        stats['alpha0'][iteration] = sample['alpha0']
        stats['gamma'][iteration] = sample['gamma']
        stats['K'][iteration] = sample['K']
        stats['stan_time'][iteration] = time_stan_finish-time_stan_start
        stats['jll'][iteration] = ihmm_joint_llike(Y = Y, rt = rt, 
                                                   correct_response = correct_response,
                                                   _filters = _filters,
                                                   Pi = sample['Pi'], Phi = sample['Phi'], 
                                                   tau0 = sample['tau0'], mu = sample['mu'], sigma = sample['sigma'], 
                                                   alpha0 = sample['alpha0'], gamma = sample['gamma'], Beta = sample['Beta'],
                                                   epsilon = epsilon,
                                                   hypers = hypers,
                                                  )

        #    if ((iteration-numb) % numi) == 0:
        print('Iteration: ', iteration, ':  JL =', stats['jll'][iteration], 'K =', sample['K'], 'epsilon =', epsilon)
        print('             tau0 =', sample['tau0'], ' mu =', sample['mu'], ' sigma =', sample['sigma'])
        print('             alpha0 =', sample['alpha0'], ' gamma =', sample['gamma'])
        print('Stan time: ', stats['stan_time'][iteration])
        if (iteration+1 >= sampling_params['numb']) :
            _sample = sample.copy()
            _sample['iteration'] = iteration
            _sample['loglike'] = stats['jll'][iteration]
            #_sample['S'] = sample['S'].copy()
            _sample['pi'] = sample['Pi'].copy()
            _sample['phi'] = sample['Phi'].copy()
            _sample['jll'] = stats['jll'][iteration]
            _sample['epsilon'] = epsilon
            
            if check_error(_sample) and len(S)>0:
                _sample = S[-1]
                print('Sample parameters are out of plausible range. Backtracking...')
                print('Sample: ', _sample)
            
            S.append(_sample.copy())

        with open(output,'wb') as file:
            pickle.dump(dict(stats = stats, 
                             samples = S, 
                             hypers = hypers,
                             sampling_params = sampling_params,
                            ),
                        file)
            file.close()
        iteration += 1

    return(S, stats)

def sample_beam_ihmm_from_file(config_filename, data_filename = None):
    config = configparser.ConfigParser()
    config.read(os.path.join(CONFIG_DIR,config_filename))

    if data_filename is not None:
        with open(os.path.join(DATA_DIR,data_filename),'rb') as file:
            data = pickle.load(file)
        output_file = os.path.join(OUTPUT_DIR,utils.get_git_hash('..')[:7]+'_'+config_filename[:-4]+'_'+data_filename[:-4]+'_samples.pkl')
        print('Importing stimulus and reaction time data from', data_filename)
        print('Writing results into', output_file)

    print('-------------------------------------------------')
    print('--- iHMM Cognitive Tomography by Balazs Torok ---')
    print('-------------------------------------------------')
    print()
    print('% Thanks to J Van Gael for providing code of iHMM Beam Sampling')
    print('% This code makes heavy use of his work')
    print()
    print('Initiating sampling using config file:', config_filename)
    
    print('# Hyperparameter settings #')
    hypers = dict()
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
    
    if 'epsilon_max' not in config['HYPERS']:
        hypers['epsilon_max'] = 0.2
    if 'epsilon_min' not in config['HYPERS']:
        hypers['epsilon_min'] = 2e-2
    
    sampling_params = dict()
    print('# Sampling parameters #')
    for key in config['SAMPLING_PARAMS']:
        sampling_params[key] = int(config['SAMPLING_PARAMS'][key])    
        print('    ',key,'=',sampling_params[key])
    print()
    
    K0 = int(config['SAMPLING_PARAMS']['K0'])
    rt = data['rt']
    for i in range(len(rt)):
        if rt[i] < 1:
            rt[i] = 1
    Y = data['Y']
    hypers['H'] = np.array([np.repeat(1,np.max(Y)+1)])*hypers['h']
    if 'correct_response' in data:
        correct_response = data['correct_response']*(np.array(rt)>1)
    else:
        correct_response = np.ones(len(rt), dtype = np.int32)*(np.array(rt)>1)
    
    _filters = filters.get_filter(data,config)
    
    sample_beam_ihmm(K0, Y, rt, correct_response, _filters, hypers, sampling_params, output = output_file)
    print('--------------------------')
    print('--- SAMPLING FINISHED. ---')
    print('--------------------------')
    print()