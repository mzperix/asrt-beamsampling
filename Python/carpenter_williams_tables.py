from asrt_table_definitions import *
from rt_model import *

import scipy.stats as sp
import numpy as np

schema = dj.schema('asrt', locals())

STAN_MODEL_FILE = 'cw_no_lapse_stan_model.pkl'

@schema
class CW_Prior(dj.Manual):
	definition = """
	prior_id: int unsigned 
	---
	tau0_alpha: double
	tau0_scale: double
	mu_alpha: double
	mu_scale: double
	sigma_alpha: double
	sigma_scale: double
	p_lapse_alpha: double
	p_lapse_beta: double
	l_alpha: double
	l_scale: double 
	"""

	def PlotPrior(key):
		rows = (CW_Prior() & key).fetch(as_dict = True)
		for row in rows:
			print('Prior id: ', row['prior_id'])
			plt.subplot(2,3,1)
			Prior.PlotGamma(row['tau0_alpha'],row['tau0_scale'])
			plt.title('tau0')
			plt.subplot(2,3,2)
			Prior.PlotGamma(row['mu_alpha'],row['mu_scale'])
			plt.title('mu')
			plt.subplot(2,3,3)
			Prior.PlotGamma(row['sigma_alpha'],row['sigma_scale'])
			plt.title('sigma')
			plt.subplot(2,3,4)
			Prior.PlotBeta(row['p_lapse_alpha'],row['p_lapse_beta'])
			plt.title('p_lapse')
			plt.subplot(2,3,5)
			Prior.PlotGamma(row['l_alpha'],row['l_scale'])
			plt.title('l')
			plt.show()
		

@schema
class CW_SyntheticParameter(dj.Manual):
	definition = """
	parameter_id: int unsigned
	---
	->Seed
	->CW_Prior
	tau0: double
	mu: double
	sigma: double
	p_lapse: double
	l: double
	"""

	def GenerateSyntheticParameter(self, key):
		last_id = np.max(CW_SyntheticParameter().fetch("parameter_id"))
		parameter_id = last_id + 1
		seed = (Seed() & key).fetch1("seed")
		prior_params = (CW_Prior() & key).fetch1(as_dict = True)

		random.seed(seed)
		tau0, mu, sigma, p_lapse, l = CarpenterWilliamsModel.sample_params_from_prior(key)

		self.insert1(dict(parameter_id = parameter_id,
						  seed = seed,
						  prior_id = prior_params['prior_id'],
						  tau0 = tau0, 
						  mu = mu, 
						  sigma = sigma,
						  p_lapse = p_lapse,
						  l = l))


@schema
class CW_HyperParameter(dj.Lookup):
	definition = """
	hyper_id: int unsigned
	---
	tau0_sigma: double
	mu_sigma: double
	sigma_sigma: double
	p_lapse_sigma: double
	l_sigma: double
	min_proportion: double # lower bounds for mu and sigma
	numi: int unsigned # number of iterations between samples
	nums: int unsigned # number of samples
	numb: int unsigned # number of burn-in iterations
	"""


@schema
class CW_RtParameterPosterior(dj.Computed):
	definition = """
	->DataFilter.Dataset
	->PredictionModel
	->CW_Prior
	->CW_HyperParameter
	->Seed
	---
	"""

	class Sample(dj.Part):
		definition = """
		->CW_RtParameterPosterior
		iteration: int unsigned # id of iteration
		---
		loglike: double # log-likelihood
		tau0: double
		mu: double
		sigma: double
		p_lapse: double
		l: double
		"""

	## METROPOLIS-HASTINGS
	def __make_tuples(self, key):
		print('Populating  seed: '+str(key['seed'])+'   filter_id: '+str(key['filter_id'])+'   dataset_id: '+str(key['dataset_id']))
		self.insert1(key)
		prior_params = (CW_Prior() & key).fetch1(as_dict = True)
		seed = (Seed() & key).fetch1(as_dict=True)['seed']
		hypers = (CW_HyperParameter() & key).fetch1(as_dict = True)
	
		for dataset in DataFilter.Dataset.GetData(key):
			random.seed(seed)
			rt = dataset['first_rt']

			logp = PredictionModel.predict(key, dataset)

	        # Initialize parameters
			loglike = -float("Inf")
			init_trial = 1
			while (loglike == -float("Inf")) & (init_trial < 100):
				model = CarpenterWilliamsModel.initialize_from_prior(key)
				loglike = model.loglike(rt, logp)
				init_trial += 1
			
			if(loglike == -float("Inf")):
				raise InitError

			for sample in model.metropolis_sampling(rt, logp, hypers, prior_params,seed):
				key['iteration'] = sample['iteration']
				key['tau0'] = sample['tau0']
				key['mu'] = sample['mu']
				key['sigma'] = sample['sigma']
				key['p_lapse'] = sample['p_lapse']
				key['l'] = sample['l']
				key['loglike'] = np.double(sample['loglike'])
				self.Sample().insert1(key)

	def _make_tuples(self, key):
		print('Populating  seed: '+str(key['seed'])+'   filter_id: '+str(key['filter_id'])+'   dataset_id: '+str(key['dataset_id']))
		self.insert1(key)
		prior_params = (CW_Prior() & key).fetch1(as_dict = True)
		seed = (Seed() & key).fetch1(as_dict=True)['seed']
		hypers = (CW_HyperParameter() & key).fetch1(as_dict = True)
		
		# import stan model
		stan_model = pickle.load(STAN_MODEL_FILE)

		for dataset in DataFilter.Dataset.GetData(key):
			rt = dataset['first_rt']
			logp = PredictionModel.predict(key, dataset)

			data = dict(N = len(rt),
						rt = rt,
						logp = logp,
						rt_max = (DataFilter() & key).fetch()['rt_max'],
						tau0_shape = 1.0,
						tau0_scale = prior_params['tau0_scale'],
						mu_shape = 1.0,
						mu_scale = prior_params['mu_scale'],
						sigma_shape = 1.0,
						sigma_scale = prior_params['sigma_scale'],
						error_sigma = 100)

			samples = stan_model.sampling(data = data, seed = seed, iter = 800, warmup = 400)
			key['iteration'] = samples.extract(['iteration'])
			key['tau0'] = samples.extract(['tau0'])
			key['mu'] = samples.extract(['mu'])
			key['sigma'] = samples.extract(['sigma'])
			key['p_lapse'] = np.repeat([0],len(key['sigma']))
			key['l'] = np.repeat([0], len(key['sigma']))
			key['loglike'] = np.double(sample['loglike'])

