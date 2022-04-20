import numpy as np
from carpenter_williams_model import CarpenterWilliamsModel
import scipy.stats as sp
import matplotlib.pyplot as plt
import samples_utils

STATE_POSTERIOR = 'STATE_POSTERIOR'
STATE_POSTERIOR_DIFFERENCE = 'STATE_POSTERIOR_DIFFERENCE'
STATE_PRIOR = 'STATE_PRIOR'
STATE_PRIOR_ENTROPIES = 'STATE_PRIOR_ENTROPIES'
STATE_POSTERIOR_ENTROPIES = 'STATE_POSTERIOR_ENTROPIES'
PRED_DIFFERENCE = 'PRED_DIFFERENCE'
PRED_PROBS = 'PRED_PROBS'
PRED_ENTROPIES = 'PRED_ENTROPIES'
LOG_PRED_PROB = 'LOG_PRED_PROB'
LOG_PRED_PROBS = 'LOG_PRED_PROBS'
PRED_PROB_ALL = 'PRED_PROB_ALL'

class IdealObserverSamples:
    """Stores and uses samples obtained from the
    posterior sampling of bsrt_sampling_rt_s_hat.
    """

    def __init__(self, **kwargs):
        if 'hypers' in kwargs:
            self.hypers = hypers
        if 'samples' in kwargs:
            self.samples = kwargs['samples']
        else:
            self.samples = samples_utils.import_samples(**kwargs)

    def instantiate_sample(self, sample):
        if 'rt_max' in sample:
            rt_max = sample['rt_max']
        else:
            rt_max = 5000
        if 'pi_markov' in sample:
            self.observer = IdealObserver(tau0 = sample['tau0'],
                                          mu = sample['mu'],
                                          sigma = sample['sigma'],
                                          Pi = sample['pi'],
                                          Phi = sample['phi'],
                                          w_markov = sample['w_markov'],
                                          pi_markov = sample['pi_markov'],
                                          rt_max = rt_max)
        else:
            self.observer = IdealObserver(tau0 = sample['tau0'],
                                          mu = sample['mu'],
                                          sigma = sample['sigma'],
                                          Pi = sample['pi'],
                                          Phi = sample['phi'],
                                          rt_max = rt_max)

    def generate_predictive_probabilities(self, Y, hypers = dict()):
        """Calculates the log_pred_prob and pred_prob_all for each sample."""
        log_pred_prob = []
        pred_prob_all = []
        for sample in self.samples:
            self.instantiate_sample(sample)
            prediction = self.observer.generate_predictive_probabilities(Y, hypers)
            log_pred_prob.append(prediction[LOG_PRED_PROB])
            pred_prob_all.append(prediction[PRED_PROB_ALL])
        return({LOG_PRED_PROBS: np.array(log_pred_prob),
                PRED_PROB_ALL: np.array(pred_prob_all)})

    def generate_marginal_predictive_probabilities(self, Y, hypers = dict()):
        """Returns a dict of log_pred_prob and pred_prob_all
        the log predicted prob for each Y and the full
        predictive distribution for each position."""
        prediction = self.generate_predictive_probabilities(Y, hypers)

        return({LOG_PRED_PROB: np.mean(np.array(prediction[LOG_PRED_PROB]), 0),
                PRED_PROB_ALL: np.mean(np.array(prediction[PRED_PROB_ALL]), 0)})

    def predict_rt(self, Y, hypers = dict()):
        rt = []
        for sample in self.samples:
            self.instantiate_sample(sample)
            prediction = self.observer.generate_predictive_probabilities(Y, hypers)
            rt.append((self.observer.tau0-prediction[LOG_PRED_PROB])/self.observer.mu)
        return(1/np.mean(1/np.array(rt),0))

    def calculate_residuals(self, data):
        """Calculates the r from the CW model
        data includes Y, correct_response and rt fields
        calculates residuals for only correct responses."""
        residuals = []
        for sample in self.samples:
            self.instantiate_sample(sample)
            prediction = self.observer.generate_predictive_probabilities(data['Y'])
            #log_pred_probs = log_pred_probs[data['correct_response']==1]
            #rt = data['rt'][data['correct_response']==1]
            rt = data['rt']
            residuals.append((self.observer.tau0-prediction[LOG_PRED_PROB])/rt)
        return(np.array(residuals))

    def seed_and_predict(self, Y):
        predictions = []
        for sample in self.samples:
            self.instantiate_sample(sample)
            predictions.append(self.observer.seed_and_predict(Y))
        rt = self.predict_rt(Y)[-1]
        return {
            'predicted_probability': predictions,
            'rt': rt}


class IdealObserver:
    def __init__(self, tau0, mu, sigma, Pi, Phi, rt_max = 5000, w_markov = 0.0, pi_markov = np.zeros((4,4))):
        self.tau0 = tau0
        self.mu = mu
        self.sigma = sigma
        self.Pi = Pi
        self.Phi = Phi
        self.rt_max = rt_max
        self.w_markov = w_markov
        self.pi_markov = pi_markov
        self.cm = CarpenterWilliamsModel(tau0 = self.tau0,
                                         mu = self.mu,
                                         sigma = self.sigma,
                                         p_lapse = 0.0,
                                         l = 500)


    def from_sample(sample):
        if 'pi_markov' in sample:
            return(IdealObserver(tau0 = sample['tau0'],
                                 mu = sample['mu'],
                                 sigma = sample['sigma'],
                                 Pi = sample['pi'],
                                 Phi = sample['phi'],
                                 rt_max = 5000))
        else:
            return(IdealObserver(tau0 = sample['tau0'],
                                 mu = sample['mu'],
                                 sigma = sample['sigma'],
                                 Pi = sample['pi'],
                                 Phi = sample['phi'],
                                 w_markov = sample['w_markov'],
                                 pi_markov = sample['pi_markov'],
                                 rt_max = 5000))

    def to_dict(self):
        return {
            'tau0': self.tau0,
            'mu': self.mu,
            'sigma': self.sigma,
            'phi': self.Phi,
            'pi': self.Pi,
            'rt_max': self.rt_max}

    def asrt_ground_truth(tau0, mu, sigma, transition_noise, emission_noise, sequence, rt_max=5000):
        assert len(sequence)==4
        sorted_sequence = np.copy(sequence)
        sorted_sequence.sort()
        assert np.array_equal(sorted_sequence, np.array([0,1,2,3]))

        Pi  = np.ones((9,9))*transition_noise
        Pi[1:,0] = 0
        Pi[0,:] = np.array([0]+[1/8]*8)
        Phi = np.ones((9,4))*emission_noise
        Phi[0,:] = np.ones(4)/4
        for s in range(1,8):
            Pi[s,s+1] += 1
            if s%2 == 0:
                Phi[s,sequence[s//2-1]] += 1
            Pi[s,:] /= np.sum(Pi[s,:])
            Phi[s,:] /= np.sum(Phi[s,:])
        Pi[-1,1] += 1
        Pi[-1,:] /= np.sum(Pi[-1,:])
        Phi[-1,sequence[-1]] += 1
        Phi[-1,:] /= np.sum(Phi[-1,:])
        return IdealObserver(tau0, mu, sigma, Pi, Phi, rt_max)

    def generate_predictive_probabilities(self, Y, hypers=dict()):
        if 'restart_stream' in hypers:
            restart_stream=hypers['restart_stream']
        else:
            restart_stream=np.zeros(len(Y))
            restart_stream[0]=1

        T = np.size(Y)
        K = np.size(self.Pi,0)
        log_pred_prob = np.zeros(T)
        pred_prob_all = np.zeros((np.size(self.Phi,1),T))
        state_posterior = np.ones((K,T))
        state_prior = np.ones((K,T))

        for t in range(0,T):
            if restart_stream[t] == 1:
                state_prior[:,t] = self.Pi[0,0:K]+1e-10
                state_prior[:,t] /= np.sum(state_prior[:,t])
            else:
                state_prior[:,t] = np.matmul(np.transpose(self.Pi[0:K,0:K]),state_posterior[:, t-1])+1e-10
                state_prior[:,t] = state_prior[:,t] / np.sum(state_prior[:,t])

            pred_prob_all[:,t] = np.matmul(self.Phi.transpose(), state_prior[:,t])
            pred_prob_all[:,t] /= np.sum(pred_prob_all[:,t])
            state_posterior[:,t] = self.Phi[:, Y[t]]*state_prior[:,t]
            state_posterior[:,t] /= np.sum(state_posterior[:,t])
            log_pred_prob[t] = np.log(pred_prob_all[Y[t],t]+1e-10)
        return {
            STATE_POSTERIOR: state_posterior,
            STATE_PRIOR: state_prior,
            LOG_PRED_PROB: log_pred_prob,
            PRED_PROB_ALL: pred_prob_all}

    def generate_rt(self,Y):
        prediction = self.generate_predictive_probabilities(Y)
        rt = self.cm.sample_rt(logp = prediction[LOG_PRED_PROB], rt_max = 5000)
        return(prediction[STATE_POSTERIOR], prediction[LOG_PRED_PROB], rt)

    def loglike(self, Y, rt):
        prediction = self.generate_predictive_probabilities(Y)
        return(self.cm.loglike(rt, prediction['LOG_PRED_PROB']))

    def seed_and_predict(self, Y):
        K = np.size(self.Pi,0)
        internal_state_distribution = np.ones(K)/K
        for y in Y[:-1]:
            internal_state_distribution *= self.Phi[:,y] # condition
            internal_state_distribution = np.matmul(
                np.transpose(self.Pi[0:K,0:K]),
                internal_state_distribution) + 1e-10 # propagate
            internal_state_distribution /= np.sum(internal_state_distribution)
        prediction = np.sum(internal_state_distribution*self.Phi[:,Y[-1]])
        return prediction


class GroundTruth(IdealObserverSamples):
    def __init__(self, transition_noise, emission_noise, sequence, **kwargs):
        super().__init__(**kwargs)
        self.transition_noise = transition_noise
        self.emission_noise = emission_noise
        self.sequence = sequence

    def instantiate_sample(self, sample):
        if 'rt_max' in sample:
            rt_max = sample['rt_max']
        else:
            rt_max = 5000
        self.observer = IdealObserver.asrt_ground_truth(
            tau0 = sample['tau0'],
            mu = sample['mu'],
            sigma = sample['sigma'],
            transition_noise = self.transition_noise,
            emission_noise = self.emission_noise,
            rt_max = rt_max,
            sequence = self.sequence)


def marginal_prediction(samples, Y):
    log_pred_probs = np.zeros((len(samples),len(Y)))
    for i, s in enumerate(samples):
        observer = IdealObserver(tau0 = s['tau0'],
                                 mu = s['mu'],
                                 sigma = s['sigma'],
                                 Pi = s['pi'],
                                 Phi = s['phi'],
                                 )
        prediction = observer.generate_predictive_probabilities(Y)
    return(np.mean(prediction[LOG_PRED_PROBS],axis=0))

def marginal_prediction_full(samples, Y):
    pred_probs = np.zeros((len(samples),np.max(Y)+1,len(Y)))
    for i, s in enumerate(samples):
        observer = IdealObserver(tau0 = s['tau0'],
                                 mu = s['mu'],
                                 sigma = s['sigma'],
                                 Pi = s['pi'],
                                 Phi = s['phi'],
                                 )
        pred_probs[i,:,:] = observer.generate_predictive_probabilities(Y)[PRED_PROB_ALL]
    return(np.mean(np.log(pred_probs+1e-10),0))

def state_and_prediction_entropy(model, Y):
    def difference(v):
        sorted_v = np.sort(v)
        return sorted_v[-1]-sorted_v[-2]

    samples = model.samples
    posterior_entropies = np.zeros((len(samples), len(Y)))
    prior_entropies = np.zeros((len(samples), len(Y)))
    state_posterior_difference = np.zeros((len(samples), len(Y)))
    pred_entropies = np.zeros(len(Y))
    pred_probs = np.zeros((len(samples), len(Y)))
    pred_prob_all = np.zeros((len(samples), np.size(samples[0]['phi'],1), len(Y)))
    pred_difference = np.zeros((len(samples), len(Y)))

    for i, s in enumerate(model.samples):
        model.instantiate_sample(s)
        d = model.observer.generate_predictive_probabilities(Y)

        # Calculate entropies each Y
        posterior_entropies[i, :] =        np.array([sp.entropy(d[STATE_POSTERIOR][:,t]) for t in range(len(Y))])
        prior_entropies[i, :] =            np.array([sp.entropy(d[STATE_PRIOR][:,t]) for t in range(len(Y))])
        state_posterior_difference[i, :] = np.array([difference(d[STATE_POSTERIOR][:,t]) for t in range(len(Y))])
        pred_probs[i, :] =                 np.array([d[PRED_PROB_ALL][Y[t],t] for t in range(len(Y))])
        pred_prob_all[i] =                 d[PRED_PROB_ALL]
        pred_difference[i, :] =            np.array([difference(d[PRED_PROB_ALL][:,t]) for t in range(len(Y))])

    pred_entropies = np.array([sp.entropy(np.mean(pred_prob_all,0)[:,t]) for t in range(len(Y))])

    return {
        STATE_POSTERIOR_ENTROPIES: posterior_entropies,
        STATE_PRIOR_ENTROPIES: prior_entropies,
        STATE_POSTERIOR_DIFFERENCE: state_posterior_difference,
        PRED_ENTROPIES: pred_entropies,
        PRED_DIFFERENCE: pred_difference,
        PRED_PROBS: pred_probs}
