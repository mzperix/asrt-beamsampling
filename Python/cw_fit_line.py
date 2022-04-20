## Line fitting function for Carpenter-Williams model

import numpy as np
import pandas as pd

from asrt_table_definitions import *
from carpenter_williams_model import *
from rt_model import *

def quantiles(x,quantiles):
	x_sorted = np.sort(x)
	result = []
	for q in quantiles:
		result.append(x_sorted[int(x_sorted[q*len(x)])])
	return(result)


def cw_fit_line(d, n_bins = 10):
	# data should contain fields: rt, logp
	
	# first let us bin logp
	data = d.copy()
    data['logp'] = PredictionModel.naive_probabilities('explicit',data.trial_type)
    data['logp_bin']  = pd.cut(data['logp'],10)
    data['logp_mean'] = data.groupby(['logp_bin'], group_keys=False).apply(lambda x: x.logp*0+x.logp.mean())
    
    for logp in data.logp_mean.unique():
    	rt = data.