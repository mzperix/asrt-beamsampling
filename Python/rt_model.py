# RT model
import numpy as np
import scipy.stats as sp
import random

import matplotlib.pyplot as plt

from asrt_table_definitions import *

schema = dj.schema('asrt', locals())

import logging
#logging.basicConfig(filename='rt_model.log',
#                    level=logging.INFO, 
#                    format='%(asctime)s %(message)s')

class Prior:
	resolution = 150

	def PlotGamma(alpha, scale, show = False):
		mean = alpha*scale
		sd = np.sqrt(alpha) * scale
		x = np.linspace(mean-3*sd, mean+8*sd, Prior.resolution)
		y = sp.gamma.pdf(x,alpha,scale = scale)
		plt.plot(x,y)
		if(show):
			plt.show()

	def PlotGauss(mu, sigma,show = False):
		x = np.linspace(mu-6*sigma,mu+6*sigma,Prior.resolution)
		y = sp.norm.pdf(x,mu,sigma)
		plt.plot(x,y)
		if(show):
			plt.show()

	def PlotBeta(alpha, beta, show = False):
		x = np.linspace(0,1,Prior.resolution)
		y = sp.beta.pdf(x,alpha, beta)
		plt.plot(x,y)
		if(show):
			plt.show()


@schema
class Seed(dj.Lookup):
	definition = """
	seed: int unsigned # value of the seed
	---
	"""


@schema
class DataFilter(dj.Manual):
	definition = """
	filter_id: int unsigned
	---
	first_acc: int unsigned # 0 if not filtering for accuracy, 1 if do
	rt_min: int unsigned # lower bound for rt in msec
	rt_max: int unsigned # upper bound for rt in msec
	"""

	class Dataset(dj.Part):
		definition = """
		->DataFilter
		dataset_id: int unsigned
		---
		->Session
		block: longblob # list of blocks to be used
		trial: longblob # list of trials to be used
		"""

		def GetData(key,cols=None):
			filter_rows = (DataFilter().Dataset() * 
				           DataFilter() & 
				           key).fetch(as_dict = True)
			for row in filter_rows:
				session_key = {"experiment_id": row['experiment_id'],
								"session_id": row['session_id'],
								"session_name": row['session_name'],
								"participant_id": row['participant_id']}
				if(cols == None):
					data = (Session().AsrtData() & session_key).fetch()
				else:
					data = (Session().AsrtData() & session_key).fetch(cols)
				data = pd.DataFrame(data)
				_data = data.loc[(data.block.isin(np.array(row['block'])) &
							 data.trial.isin(np.array(row['trial'])) &
							 (data.first_acc >= row['first_acc']) &
							 (data.first_rt >= row['rt_min']) &
							 (data.first_rt <= row['rt_max']))]
				yield _data


@schema
class PredictionModel(dj.Lookup):
	definition = """
	prediction_model_id: int unsigned
	---
	prediction_model_name: enum("naive","triplet","ihmm") # name of the prediction model
	implicit_explicit: enum("implicit","explicit") # which asrt type it corresponds to
	"""

	def naive_probabilities(imp_exp, trial_type):

		if (imp_exp == 'implicit'):
			log25 = np.log(0.25)
			logp = np.repeat(log25, len(trial_type))

		if (imp_exp == 'explicit'):
			log1 = 0
			log25 = np.log(0.25)
			logp = np.zeros(len(trial_type))
			for i,trial in enumerate(trial_type):
				if (trial == 'P'):
					logp[i] = log1
				if (trial == 'R'):
					logp[i] = log25
				if (trial == 'Prac'):
					logp[i] = log25
		return(logp)

	def predict(key, dataset):
		prediction_model = (PredictionModel() & key).fetch1(as_dict = True)

		if ((prediction_model['prediction_model_name'] == 'naive') & 
			(prediction_model['implicit_explicit'] == 'explicit')):
			logp = PredictionModel.naive_probabilities('explicit', dataset['trial_type'])
			return(logp)

		raise NotImplementedError('This prediction model is not implemented.')


