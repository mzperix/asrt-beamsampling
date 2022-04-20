import datajoint as dj

from table_definitions import *
import rt_model

schema = dj.schema('asrt', locals())

def filtered_fetch(key):
	return(((Session.AsrtData() & key)*(AsrtFilter() & key) & 
		'first_acc >= acc' &
		'first_rt >= min_rt' &
		'trial > start_random * 5').fetch())


@schema
class AsrtFilter(dj.Lookup):
	definition = """
	filter_id: tinyint unsigned
	---
	min_rt: int   # lower bound for rt
	max_rt: int 

	acc: tinyint # 1 if filter for accurate trials
	start_random: tinyint # 1 if filter for randoms at the start
	"""

def Filter_session(session, filt):
	pass


@schema
class RtMcmcParam(dj.Lookup):
	definition = """
	# MCMC parameter settings for rt_model
	parameter_id: int
	seed: int
	---
	tau0: float # Initial parameters
	mu: float
	sigma: float
	p_lapse: float
	l: float
	
	tau0_sigma: float
	mu_sigma: float
	sigma_sigma: float
	p_lapse_sigma: float
	l_sigma: float

	maxiter: int
	"""


@schema
class RtParamSamples(dj.Lookup):
	definition = """
	# MCMC samples for rt_model parameters
	->Session
	->RtMcmcParam
	->AsrtFilter
	---
	"""

	class Stats(dj.Part):
		definition = """
		->RtParamSamples
		iteration: int
		---
		tau0: float
		mu: float
		sigma: float
		p_lapse: float
		l: float
		loglike: float
		"""

	def _make_tuples(self, key):
		rt = (asrt.Session() & key).fetch('first_rt')