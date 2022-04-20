from asrt_table_definitions import *
from rt_model import *
from carpenter_williams_model import *

import numpy as np
import plots
import matplotlib.pyplot as plt

from IPython.display import display

def sorted_groupby(df, groupby):
    _df = df.sort(groupby)
    start = 0
    prev = _df[groupby].iloc[start]
    for i, x in enumerate(_df[groupby]):
        if x != prev:
            yield prev, _df.iloc[start:i]
            prev = x
            start = i
    # need to send back the last group
    yield prev, _df.iloc[start:]

class Plots:
	pass

class CW_PosteriorPlots(Plots):

	def plot_posterior_samples(key):

		samples = pd.DataFrame((CW_RtParameterPosterior.Sample() & key).fetch())

		for (filter_id, dataset_id, prediction_model_id, prior_id, hyper_id), df in samples.groupby(['filter_id', 'dataset_id', 'prediction_model_id', 'prior_id', 'hyper_id']):

			for dataset in DataFilter.Dataset.GetData(dict(filter_id = filter_id, dataset_id = dataset_id)):
				rt = dataset.first_rt
				participant_id = dataset.iloc[0].participant_id
				plt.plot(rt, '.')
				plt.title('First rt,   Participant: '+str(participant_id))
				plt.show()

				hypers = (CW_HyperParameter() & dict(hyper_id = hyper_id)).fetch(limit = 1)
				display(pd.DataFrame(hypers))
				
				sorted_df = df.sort_values(by = 'iteration')
				
				n_seed_shown = len(df.seed.unique())
				table = df[df.iteration == df.iteration.max()]
				display(table)
				

				for i, seed in enumerate(df.seed.unique()):
					sample = sorted_df.iloc[-n_seed_shown+i]
					#print('Hypers: ', hypers)
					title = ('Participant: '+str(participant_id)+
							'   loglike: '+str(sample.loglike)+'   iteration: '+str(sample.iteration))
					model = CarpenterWilliamsModel(**sample.to_dict())
					logp = PredictionModel.predict(dict(prediction_model_id = prediction_model_id), dataset)
					rt_sample = model.sample_rt(np.tile(logp,30))

					plt.subplot(n_seed_shown,2,1+2*i)
					_, bins, _ = plt.hist([rt[rt<1000],rt_sample[rt_sample < 1000]], 35, normed = True)
					#plt.hist(rt_sample, bins = bins, normed = True)
					plt.title(title)
					
					plt.subplot(n_seed_shown,2,2+2*i)
					plots.QQ_plot(rt, rt_sample)
					plt.show()


				for seed in samples.seed.unique():
					x = range(len(df[df.seed == seed].seed))
					plt.subplot(3,3,1)
					plt.plot(x,df[df.seed == seed].tau0)
					plt.title('tau0')

					plt.subplot(3,3,2)
					plt.plot(x,df[df.seed == seed].mu)
					plt.title('mu')

					plt.subplot(3,3,3)
					plt.plot(x,df[df.seed == seed].sigma)
					plt.title('sigma')

					plt.subplot(3,3,4)
					plt.plot(x,df[df.seed == seed].p_lapse)
					plt.title('p_lapse')

					plt.subplot(3,3,5)
					plt.plot(x,df[df.seed == seed].l)
					plt.title('l')

					plt.subplot(3,3,6)
					plt.plot(x,df[df.seed == seed].loglike)
					plt.title('loglike')

					plt.subplot(3,3,9)
					plt.plot(x[200:],df[df.seed == seed].loglike[200:])
					plt.title('loglike')
				plt.show()