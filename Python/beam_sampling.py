### BEAM SAMPLING FUNCTIONS

### Based on: Matlab beamsampling functions

import numpy as np
import scipy.stats as sp
from scipy.special import gammaln

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

# Matlab function: dirichlet_sample.m
def sample_dirichlet(alpha):
    sample = random_gamma(alpha, 1.0)
    sample = sample / np.sum(sample)
    return(sample)


# Matlab function: SampleEmissionMatrix(S, Y, K, H)
def sample_emission_matrix(S, Y, K, H):
    # S: state sequence (np.array)
    # Y: emission sequence (np.array)
    # K: number of represented states
    # H: parameter of emission prior (~ Dir(H))

    L = np.size(H, 1)
    Phi = np.zeros((K,L))

    for k in range(K):
        empirical = np.zeros((1,L))
        for l in range(L):
            empirical[0,l] = np.sum(Y[S==k]==l)
        Phi[k,:] = sample_dirichlet(empirical+H)
    return(Phi)


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
def ihmm_joint_llike(S, Y, Beta, alpha0, H):
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
def sample_beam_ihmm(Y, hypers, numb, nums, numi, S0):
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

    # Sample the emission and transition probabilities.
    sample['Phi'] = sample_emission_matrix(sample['S'], Y, sample['K'], hypers['H'])
    sample['Pi'] = sample_transition_matrix(sample['S'], sample['alpha0']*sample['Beta'])
    sample['Pi'] = np.delete(sample['Pi'],sample['K'],0)

    iteration = 0;
    # print values
    
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
        dyn_prog = np.zeros((sample['K'],T))

        dyn_prog[:, 0] = sample['Pi'][0,0:sample['K']] > u[0]
        stats['trellis'][iteration] += np.sum(np.sum(dyn_prog[0,:]))

        for k in range(sample['K']):
            dyn_prog[k,0] = sample['Phi'][k, Y[0]] * dyn_prog[k,0]
        dyn_prog[:,0] = dyn_prog[:,0] / np.sum(dyn_prog[:,0])

        for t in range(1,T):
            A = sample['Pi'][0:sample['K'], 0:sample['K']] > u[t]
            dyn_prog[:,t] = np.matmul(np.transpose(A),dyn_prog[:,t-1])
            stats['trellis'][iteration] += np.sum(np.sum(A))
            assert(np.size(dyn_prog,0) == sample['K'])
            assert(np.size(sample['Phi'], 0) == sample['K'])
            for k in range(sample['K']):
                dyn_prog[k,t] = sample['Phi'][k, Y[t]] * dyn_prog[k, t]
            dyn_prog[:,t] = dyn_prog[:,t] / np.sum(dyn_prog[:,t])
            
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

            # Resample Beta given the transition probabilities.
            sample['Beta'], sample['alpha0'], sample['gamma'], _, _ = sample_ihmm_hyper(sample['S'], sample['Beta'], sample['alpha0'], sample['gamma'], hypers, 20)

            # Resample the Phi's given the new state sequence
            sample['Phi'] = sample_emission_matrix(sample['S'], Y, sample['K'], hypers['H'])

            # Resample transition probabilities.
            sample['Pi'] = sample_transition_matrix(sample['S'], sample['alpha0'] * sample['Beta'])
            sample['Pi'] = np.delete(sample['Pi'], sample['K'], axis = 0)

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
            stats['jll'][iteration] = ihmm_joint_llike(sample['S'], Y, sample['Beta'], sample['alpha0'], hypers['H'])

            if (iteration % 200) == 0:
                print('Iteration: ', iteration, ':  K =', sample['K'], ' alpha0 =', sample['alpha0'],
                      ' gamma =', sample['gamma'], ' JL =', stats['jll'][iteration])
            if (iteration+1 >= numb) and (((iteration+1-numb) % numi) == 0):
                S.append(sample.copy())

            iteration += 1
        else:
            print('Wasted computation as there were no paths through the iHMM.')

    return(S, stats)