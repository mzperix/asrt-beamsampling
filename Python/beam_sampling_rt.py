### BEAM SAMPLING FUNCTIONS

### CONDITIONING ON RT !!!!


### Based on: Matlab beamsampling functions

import numpy as np
import scipy.stats as sp
from scipy.special import gammaln
import carpenter_williams_model as cwm
import pickle
import pystan
import time
import utils
import sys, os

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
sm_emission_matrix = pickle.load(open(utils.find_file('sm_emission_matrix.pkl'),'rb'))
sm_rt_model = pickle.load(open(utils.find_file('sm_rt_model.pkl'),'rb'))
sm_joint_e_rt = pickle.load(open(utils.find_file('sm_joint_e_rt.pkl'),'rb'))


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
def ihmm_joint_llike(S, Y, Beta, alpha0, H, rt_params, rt, Phi):
    # S: state sequence (np.array)
    # Y: emission sequence (np.array)
    # Beta: parameter of transition matrix prior (result of DP)
    # alpha0: hyperparameter
    # H: parameter of dirichlet emission prior

    K = np.max(S)+1
    L = np.size(H,1)
    T = np.size(S)
    N = np.zeros((K,K))
    E = np.zeros((K,L))

    # Compute all transition and emission counts.
    N[0,S[0]] = 1
    E[0,Y[0]] += 1
    for t in range(1,T):
        N[S[t-1],S[t]] += 1
        E[S[t], Y[t]] += 1

    # Compute the log likelihood
    logp = 0
    rt_model = cwm.CarpenterWilliamsModel(tau0 = rt_params['tau0'],
                                         mu = rt_params['mu'],
                                         sigma = rt_params['sigma'],
                                         p_lapse = 0.0,
                                         l = 1)

    log_phi = np.log(1e-100+Phi)
    log_pred_prob = np.zeros(T)
    for t in range(T):
        log_pred_prob[t] = log_phi[S[t],Y[t]]

    for k in range(K):
        R = np.append(N[k,:],0)+alpha0*Beta;
        ab = alpha0*Beta;
        nzind = (R != 0)
        # Add transition likelihood
        logp = logp + gammaln(alpha0)
        logp = logp - np.sum(gammaln(np.sum(np.append(N[k,:],0))+alpha0))
        logp = logp + np.sum(gammaln(R[nzind]))
        logp = logp - np.sum(gammaln(ab[nzind]))

        # Add emission likelihood
        logp = logp + (gammaln(np.sum(H))
                    - np.sum(gammaln(H))
                    + np.sum(gammaln(H+E[k,:]))
                    - gammaln(np.sum(H+E[k,:])))

        # Add rt likelihood
        logp = logp + np.sum(rt_model.loglike(rt, log_pred_prob))

    return(logp)


# Matlab function: iHmmHyperSample
# Resamples hyperparameters of an infinite hmm
def sample_ihmm_hyper(S, ibeta, ialpha0, igamma, hypers, numi):
#   [sbeta, salpha0, sgamma, N, M] = ...
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
    
    K = np.size(ibeta,1)-1 # number of states in iHmm
    T = np.size(S)       # length of iHmm sequence

    # Computer N: state transition counts
    N = np.zeros((K,K), dtype=int)
    N[0,S[0]] = 1
    for t in range(1,T):
        N[S[t-1],S[t]] += 1

    # Compute M: number of similar dishes in each restaurant
    M = np.zeros((K,K),dtype = int)
    for j in range(K):
        for k in range(K):
            if N[j,k] == 0:
                M[j,k] = 0
            else:
                for l in range(N[j,k]):
                    M[j,k] += (random_uniform() < ((ialpha0 * ibeta[0,k]) / (ialpha0 * ibeta[0,k] + l)))

    # Resample beta
    ibeta = sample_dirichlet([np.append(np.sum(M,0), igamma)])

    # Resample alpha
    if 'alpha0' in hypers:
        ialpha0 = hypers['alpha0']
    else:
        for iteration in range(numi):
            try: 
                w = random_beta(ialpha0+1, np.sum(N,1))
            except:
                print('setdiff:',np.setdiff1d(range(K),np.unique(S)))
            p = np.sum(N,1) / ialpha0
            p = p / (p+1)
            s = random_binomial(1.0, p)
            ialpha0 = random_gamma(hypers['alpha0_a'] + np.sum(np.sum(M)) - np.sum(s), 1.0 / (hypers['alpha0_b']-np.sum(np.log(w))))

    # Resample gamma (using Escobar & West, 1995)
    if 'gamma' in hypers:
        igamma = hypers['gamma']
    else:
        k = np.size(ibeta)
        m = np.sum(np.sum(M))
        for iteration in range(numi):
            mu = random_beta(igamma+1, m)
            pi_mu = 1 / (1 + (m * (hypers['gamma_b']-np.log(mu))) / (hypers['gamma_a'] + k - 1))
            if random_uniform() < pi_mu:
                igamma = random_gamma(hypers['gamma_a'] + k, 1.0 / (hypers['gamma_b'] - np.log(mu)))
            else:
                igamma = random_gamma(hypers['gamma_a'] + k - 1, 1.0 / (hypers['gamma_b'] - np.log(mu)))

    sbeta = ibeta
    salpha0 = ialpha0
    sgamma = igamma

    return(sbeta, salpha0, sgamma, N, M)


# Matlab function: iHmmSampleBeam(Y, hypers, numb, nums, numi, S0)
def sample_beam_ihmm(S0, Y, rt, hypers, numb, nums, numi, numi_rt, numi_phi):
# % IHMMSAMPLEBEAM Samples states from the iHMM with multinomial output
# % using the Beam sampler.
# %
# % [S, stats] = iHmmSampleBeam(Y, hypers, numb, nums, numi, S0) uses the
# % beam sampling training algorithm for the infinite HMM.
# %
# %   Input Parameters:
# %   - Y: training sequence of arbitrary length,
# %   - hypers: a structure that describes the hyperparameters for the beam
# %             sampler. If this structure contains alpha0 and gamma, it will
# %             not resample these during sampling. If these are not
# %             specified, one needs to specify hyperparameters for alpha0
# %             and gamma (alpha0_a, alpha0_b, gamma_a, gamma_b). hypers
# %             should also contain a prior for the emission alphabet in the
# %             field H,
# %   - numb: the number of burnin iterations,
# %   - nums: the number of samples to output,
# %   - numi: the number of sampling, iterations between two samples,
# %   - S0: is the initial assignment to the sequence.
# %
# %   Output Parameters:
# %   - S: is a cell array of sample structures where each sample contains the
# %        hidden state sequence S, the number of states K, the Beta, Pi,
# %        Phi's used for that sample.
# %   - stats: is a structure that contains a variety of statistics for every
# %            iteration of the sampler: K, alpha0, gamma, the size of the
# %            trellis and the marginal likelihood.
    
    # Initialize the sampler.
    T = np.size(Y)

    sample = dict(S = S0, 
                  K = np.max(S0)+1)

    # Setup structures to store the output
    S = []
    stats = dict(K = np.zeros(numb + (nums-1)*numi+1),
                 alpha0 = np.zeros(numb + (nums-1)*numi+1),
                 gamma  = np.zeros(numb + (nums-1)*numi+1),
                 jll = np.zeros(numb + (nums-1)*numi+1),
                 trellis = np.zeros(numb + (nums-1)*numi+1))

    # Initialize hypers; resample a few times as our initial guess might be off.
    if 'alpha0' in hypers:
        sample['alpha0'] = hypers['alpha0']
    else:
        sample['alpha0'] = random_gamma(hypers['alpha0_a'], 1.0 / hypers['alpha0_b'])

    if 'gamma' in hypers:
        sample['gamma'] = hypers['gamma']
    else:
        sample['gamma'] = random_gamma(hypers['gamma_a'], 1.0 / hypers['gamma_b'])

    for i in range(5):
        sample['Beta'] = np.ones((1, sample['K']+1)) / (sample['K']+1)
        sample['Beta'], sample['alpha0'], sample['gamma'], _, _ = sample_ihmm_hyper(sample['S'], sample['Beta'], sample['alpha0'], sample['gamma'], hypers, 20)

    # Initialize rt parameters
    sample['tau0'] = random_gamma(hypers['tau0_a'], 1.0 / hypers['tau0_b'])
    sample['mu'] = np.mean((sample['tau0']-np.log(0.25)) / rt)
    sample['sigma'] = 4*np.std(((sample['tau0']-np.log(0.25)) / rt))
    print('RT params initialized as:  tau0:', sample['tau0'], ' mu:', sample['mu'], ' sigma:', sample['sigma'])
    #sample['mu'] = random_gamma(hypers['mu_a'], 1.0 / hypers['mu_b'])
    #sample['sigma'] = random_gamma(hypers['sigma_a'], 1.0 / hypers['sigma_b'])


    # Sample the emission and transition probabilities.
    phi0 = np.ones((sample['K'],4)) / 4 # CORRECT IT LATER
    rt_params = dict(tau0 = sample['tau0'], mu = sample['mu'], sigma = sample['sigma'])
    sample['Phi'] = sample_emission_matrix(sample['S'], Y, sample['K'], hypers['H'], rt, rt_params.copy(), phi0, numi_phi)
    sample['Pi'] = sample_transition_matrix(sample['S'], sample['alpha0']*sample['Beta'])
    sample['Pi'] = np.delete(sample['Pi'],sample['K'],0)


    iteration = 0;
    # print values
    wasted_computations = 0

    while iteration <= (numb + (nums-1)*numi):
        # Safety check
        assert(np.size(sample['Phi'],0) == np.size(sample['Beta'],1) - 1)

        # Reset the trellis size count in case the previous iteration did not
        # return a samplable path.
        stats['trellis'][iteration] = 0

        # Sample the auxilary variables and extend Pi and Phi if necessary.
        u = np.zeros(T)
        for t in range(T):
            if t == 0:
                u[t] = random_uniform() * sample['Pi'][0, sample['S'][t]]
            else:
                u[t] = random_uniform() * sample['Pi'][sample['S'][t-1], sample['S'][t]]

        while np.max(sample['Pi'][:, -1]) > np.min(u): # Break the Pi[k] stick some more.
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
                #print(u)
            else:
                pg = random_beta(a, b)
            pg = np.array([np.float64(x) for x in pg])
            sample['Pi'][:, pl-1] = pg * pe
            sample['Pi'] = np.hstack([sample['Pi'], np.transpose(np.array([(1.0-pg) * pe]))])

        sample['K'] = np.size(sample['Pi'],0)

        # Safety check
        assert(sample['K'] == np.size(sample['Beta'])-1)
        assert(sample['K'] == np.size(sample['Phi'],0))

        # Resample the hidden state sequence.

        rt_model = cwm.CarpenterWilliamsModel(tau0 = sample['tau0'], mu = sample['mu'], sigma = sample['sigma'], p_lapse = 0.0, l = 1)

        start_dyn_prog = time.time()

        dyn_prog = np.zeros((sample['K'],T))

        dyn_prog[:, 0] = sample['Pi'][0,0:sample['K']] > u[0]
        stats['trellis'][iteration] += np.sum(np.sum(dyn_prog[0,:]))

        log_phi = np.log(1e-40+sample['Phi'])

        for k in range(sample['K']):
            dyn_prog[k,0] = np.exp(rt_model.loglike(rt[0],log_phi[k, Y[0]])) * sample['Phi'][k, Y[0]] * dyn_prog[k,0] 
        print(np.sum(dyn_prog[:,0]))
        dyn_prog[:,0] = dyn_prog[:,0] / np.sum(dyn_prog[:,0])

        # Precompute rt emission probabilities
        rt_probs = np.zeros((sample['K'],np.size(Y)))
        for k in range(sample['K']):
            rt_probs[k,:] = np.exp(rt_model.loglike(rt, log_phi[k, Y]))

        for t in range(1,T):
            A = sample['Pi'][0:sample['K'], 0:sample['K']] > u[t]
            dyn_prog[:,t] = np.matmul(np.transpose(A),dyn_prog[:,t-1])
            stats['trellis'][iteration] += np.sum(np.sum(A))
            assert(np.size(dyn_prog,0) == sample['K'])
            assert(np.size(sample['Phi'], 0) == sample['K'])
            for k in range(sample['K']):
                dyn_prog[k,t] = rt_probs[k,t] * sample['Phi'][k, Y[t]] * dyn_prog[k, t]
            dyn_prog[:,t] = dyn_prog[:,t] / np.sum(dyn_prog[:,t])
        
        #print('dyn prog part time: ', time.time()-start_dyn_prog)

        start_state_sample = time.time()
        # Backtrack to sample a path through the HMM.
        if (np.sum(dyn_prog[:,T-1]) != 0.0) and np.isfinite(np.sum(dyn_prog[:,T-1])):
            sample['S'][T-1] = np.sum(random_uniform() > np.cumsum(dyn_prog[:,T-1]))
            for t in range(T-2,-1,-1):
                r = dyn_prog[:,t] * (sample['Pi'][:, sample['S'][t+1]] > u[t+1])
                r = r / np.sum(r)
                sample['S'][t] = np.sum(random_uniform() > np.cumsum(r))
            # Safety check.
            assert(not np.isnan(np.sum(sample['S'][t])))

            # Cleanup our state space by removing redundant states.
            zind = np.sort(np.setdiff1d(range(sample['K']),np.unique(sample['S'])))
            zind = zind[::-1]
            for i in zind:
                sample['Beta'][0,-1] += sample['Beta'][0,i]
                sample['Beta'] = np.delete(sample['Beta'], i, axis = 1)
                sample['Pi'] = np.delete(sample['Pi'], i, axis = 0)
                sample['Pi'] = np.delete(sample['Pi'], i, axis = 1)
                sample['Phi'] = np.delete(sample['Phi'], i, axis = 0)
                sample['S'][sample['S']>i] -= 1
            sample['K'] = np.size(sample['Pi'],0)

            #print('state sampling: ', time.time()-start_state_sample)

            # Resample Beta given the transition probabilities.
            sample['Beta'], sample['alpha0'], sample['gamma'], _, _ = sample_ihmm_hyper(sample['S'], sample['Beta'], sample['alpha0'], sample['gamma'], hypers, 20)

            
            # Resample transition probabilities.
            sample['Pi'] = sample_transition_matrix(sample['S'], sample['alpha0'] * sample['Beta'])
            sample['Pi'] = np.delete(sample['Pi'], sample['K'], axis = 0)

            # Resample rt parameters AND Phi
            rt_params = dict(tau0 = sample['tau0'], mu = sample['mu'], sigma = sample['sigma'])
            #print('loglike before sampling Phi and rt params: ', ihmm_joint_llike(sample['S'], Y, sample['Beta'], sample['alpha0'], hypers['H'], rt_params, rt, sample['Phi']))
            sample['tau0'], sample['mu'], sample['sigma'], sample['Phi'] = sample_joint_e_rt(sample['S'], Y, sample['K'], hypers['H'], 
                                                                                             rt = rt, hypers = hypers, 
                                                                                             phi0 = sample['Phi'], rt_params = rt_params, 
                                                                                             numi = numi_rt)
            rt_params = dict(tau0 = sample['tau0'], mu = sample['mu'], sigma = sample['sigma'])
            #print('loglike after  sampling Phi and rt params: ', ihmm_joint_llike(sample['S'], Y, sample['Beta'], sample['alpha0'], hypers['H'], rt_params, rt, sample['Phi']))
            
            # Safety checks.
            assert(np.size(sample['Pi'],0) == sample['K'])
            assert(np.size(sample['Pi'],1) == sample['K']+1)
            assert(sample['K'] == np.size(sample['Beta'])-1)
            assert(np.min(np.min(sample['Pi'] >= 0)))
            assert(sample['K'] == np.max(sample['S'])+1)

            # Prepare next iteration.
            stats['alpha0'][iteration] = sample['alpha0']
            stats['gamma'][iteration] = sample['gamma']
            stats['K'][iteration] = sample['K']
            stats['jll'][iteration] = ihmm_joint_llike(sample['S'], Y, sample['Beta'], sample['alpha0'], hypers['H'], rt_params, rt, sample['Phi'])

            if ((iteration-numb) % numi) == 0:
                print('Iteration: ', iteration, ':  K =', sample['K'], ' alpha0 =', sample['alpha0'],
                      ' gamma =', sample['gamma'], ' JL =', stats['jll'][iteration], ' tau0 =', sample['tau0'])
            if (iteration+1 >= numb) and (((iteration+1-numb) % numi) == 0):
                _sample = sample.copy()
                _sample['iteration'] = iteration
                _sample['loglike'] = stats['jll'][iteration]
                _sample['S'] = sample['S'].copy()
                S.append(_sample.copy())

            iteration += 1
            wasted_computations = 0

        else:
            wasted_computations += 1
            print('Wasted computation as there were no paths through the iHMM.')
            if (wasted_computations > 100):
                break

    return(S, stats)