#from rt_model import *

import scipy.stats as sp
import numpy as np

STAN_MODEL_FILE = 'cw_no_lapse_stan_model.pkl'

class CarpenterWilliamsModel:

	r_min = 0.000001

	def __init__(self, **kwargs):
		self.mu = kwargs['mu']
		self.sigma = kwargs['sigma']
		self.tau0 = kwargs['tau0']
		self.l = kwargs['l']
		self.p_lapse = kwargs['p_lapse']
		#self.mean_ndt = kwargs['mean_ndt'] # mean of non-decision time
		#self.width_ndt = kwargs['width_ndt'] # width of non-decision time
		#self.rt_min = kwargs['rt_min']
		#self.rt_max = kwargs['rt_min']
		self.stats = []

	def _loglike(rt, logp, tau0, mu, sigma, p_lapse, l):
		assert np.size(rt) == np.size(logp)
		llike = 1 / rt**2 * (1-p_lapse) * sp.norm.pdf(1 / rt, mu / (tau0-logp),sigma / (tau0-logp))
		llike += np.clip(p_lapse * (sp.norm.cdf(1 / (rt - l), mu / (tau0-logp), sigma / (tau0-logp))-
    								  sp.norm.cdf(1 / rt, mu / (tau0-logp), sigma / (tau0-logp))) / l, 0, None)
		return(np.sum(np.log(llike)))

	def __loglike(rt, logp, tau0, mu, sigma, p_lapse, l):
		assert np.size(rt) == np.size(logp)
		llike1 = -2*np.log(rt) + np.log(1-p_lapse) + (sp.norm.logpdf(1 / rt, mu / (tau0-logp),sigma / (tau0-logp)) - 
		        sp.norm.logcdf(0,-mu / (tau0-logp), sigma / (tau0-logp)))
		if (p_lapse == 0.0):
			return(llike1)

		llike2 = np.clip(p_lapse * (sp.norm.cdf(1 / (rt - l), mu / (tau0-logp), sigma / (tau0-logp))-
    								  sp.norm.cdf(1 / rt, mu / (tau0-logp), sigma / (tau0-logp))) / l, 0, None)

		return(np.sum(np.log(np.exp(llike1)+llike2)))


	def loglike(self, rt, logp):
		return(CarpenterWilliamsModel.__loglike(rt, logp, self.tau0, self.mu, self.sigma, self.p_lapse, self.l))


	def logprior(prior_params, tau0, mu, sigma, p_lapse, l):
		lprior = (sp.gamma.logpdf(tau0, prior_params['tau0_alpha'], scale = prior_params['tau0_scale']) +
				  sp.gamma.logpdf(mu, prior_params['mu_alpha'], scale = prior_params['mu_scale']) +
				  sp.gamma.logpdf(sigma, prior_params['sigma_alpha'], scale = prior_params['sigma_scale']) +
				  sp.beta.logpdf(p_lapse, prior_params['p_lapse_alpha'], prior_params['p_lapse_beta']) +
				  sp.gamma.logpdf(l, prior_params['l_alpha'], scale = prior_params['l_scale']))
		return(lprior)


	def sample_params_from_prior(key):
		prior_params = (CW_Prior() & key).fetch1(as_dict = True)
		tau0 = np.double(sp.gamma.rvs(prior_params['tau0_alpha'],scale = prior_params['tau0_scale'],size = 1))
		mu = np.double(sp.gamma.rvs(prior_params['mu_alpha'],scale = prior_params['mu_scale'],size = 1))
		sigma = np.double(sp.gamma.rvs(prior_params['sigma_alpha'],scale = prior_params['sigma_scale'],size = 1))
		p_lapse = np.double(sp.beta.rvs(prior_params['p_lapse_alpha'], prior_params['p_lapse_beta'],size =1))
		l = np.double(sp.gamma.rvs(prior_params['l_alpha'],scale = prior_params['l_scale'],size =1))
		return(tau0, mu, sigma, p_lapse, l)


	def initialize_from_prior(key):
		tau0, mu, sigma, p_lapse, l = CarpenterWilliamsModel.sample_params_from_prior(key)
		return(CarpenterWilliamsModel(tau0 = tau0, mu = mu, sigma = sigma, p_lapse = p_lapse, l = l))


	def sample_rt(self, logp, rt_max = 0, plot = False):
		N = np.size(logp)
		if(rt_max > 0):
			r = np.array([sp.truncnorm.rvs(a=((self.tau0-lp)/rt_max-self.mu)/self.sigma, b=1000/self.sigma) * self.sigma + self.mu for lp in logp])
		else:
			r = sp.truncnorm.rvs(a=-self.mu/self.sigma, b=1000/self.sigma, size = N) * self.sigma + self.mu
		z = np.random.binomial(1,self.p_lapse, N)
		u = np.random.uniform(0,self.l,N)
		#ndt = np.random.uniform(self.mean_ndt-self.width_ndt / 2, self.mean_ndt + self.width_ndt / 2)
		rt = (self.tau0-logp) / r + u * z # + ndt

		if(plot):
			import matplotlib.pyplot as plt
			plt.hist(rt,50, normed = True)
			x = np.linspace(np.min(rt), np.max(rt),200)
			p = [np.exp(self.loglike(x[i],logp[i])) for i in range(len(x))]
			plt.plot(x,p)
			plt.show()
		return(rt)

	def propose_rt_params(self, hypers):
		min_proportion = hypers['min_proportion']
		_mu = -1
		_sigma = -1
		_p_lapse = -1
		_l = -1
		while (_mu < hypers['mu_sigma']*min_proportion):
			_mu = sp.norm.rvs(self.mu, hypers['mu_sigma'],1)
		while (_sigma < hypers['sigma_sigma']*min_proportion):
			_sigma = sp.norm.rvs(self.sigma, hypers['sigma_sigma'],1)
		while ((_p_lapse < 0) or (_p_lapse > 1)):
			_p_lapse = sp.norm.rvs(self.p_lapse, hypers['p_lapse_sigma'],1)
		while (_l < 0):
			_l = sp.norm.rvs(self.l, hypers['l_sigma'],1)

		_tau0 = sp.norm.rvs(self.tau0 * _mu / self.mu, hypers['tau0_sigma'],1)

		q = (sp.norm.logcdf(-hypers['mu_sigma']*min_proportion, -self.mu, hypers['mu_sigma']) - 
			 sp.norm.logcdf(-hypers['mu_sigma']*min_proportion, -_mu, hypers['mu_sigma']))
		q += (sp.norm.logcdf(-hypers['sigma_sigma']*min_proportion, -self.sigma, hypers['sigma_sigma']) -
			  sp.norm.logcdf(-hypers['sigma_sigma']*min_proportion, -_sigma, hypers['sigma_sigma']))
		q += (np.log(sp.norm.cdf(1, self.p_lapse, hypers['p_lapse_sigma']) -
			         sp.norm.cdf(0, self.p_lapse, hypers['p_lapse_sigma'])))
		q -= (np.log(sp.norm.cdf(1, _p_lapse, hypers['p_lapse_sigma'])- 
					 sp.norm.cdf(0, _p_lapse, hypers['p_lapse_sigma'])))
		q += (sp.norm.logcdf(0, -self.l, hypers['l_sigma']) -
			  sp.norm.logcdf(0, -_l, hypers['l_sigma']))

		return(np.double(_tau0), np.double(_mu), np.double(_sigma), np.double(_p_lapse), np.double(_l), q)


	def metropolis_step(self, rt, logp, loglike, hypers, prior_params):
		_tau0, _mu, _sigma, _p_lapse, _l, q = self.propose_rt_params(hypers)
		_loglike = CarpenterWilliamsModel._loglike(rt, logp, _tau0, _mu, _sigma, _p_lapse, _l)


		if (sp.uniform.rvs() > np.exp(loglike - _loglike - q + 
			CarpenterWilliamsModel.logprior(prior_params, self.tau0, self.mu, self.sigma, self.p_lapse, self.l)-
		    CarpenterWilliamsModel.logprior(prior_params, _tau0, _mu, _sigma, _p_lapse, _l))):
			self.tau0 = _tau0
			self.mu = _mu
			self.sigma = _sigma
			self.p_lapse = _p_lapse
			self.l = _l 
			return(_loglike)
		else:
			return(loglike)


	def metropolis_sampling(self, rt, logp, hypers, prior_params, seed):
		iteration = 1
		_tau0 = self.tau0
		_mu = self.mu
		_sigma = self.sigma
		_p_lapse = self.p_lapse
		_l = self.l
		_loglike = float("-inf")

		numb = hypers['numb']
		numi = hypers['numi']
		nums = hypers['nums']

		np.random.seed(seed)
		while iteration < numb + numi * nums:
			_loglike = self.metropolis_step(rt, logp, _loglike, hypers, prior_params)
			if (iteration > numb) & ((iteration - numb) % (numi) == 0):
				if (iteration - numb) % (numi*20) == 0:
					print('Iteration:  ', iteration, '   loglike:  ', _loglike)
					logging.info('Metropolis-Hastings iteration: '+str(iteration)+'  loglike: '+str(_loglike))
				
				self.stats.append(dict(hypers,
				              iteration = iteration,
				              tau0 = self.tau0,
						      mu = self.mu,
						      sigma = self.sigma,
						      p_lapse = self.p_lapse,
					    	  l = self.l,
					    	  loglike = _loglike))
				yield self.stats[-1]
			iteration += 1

