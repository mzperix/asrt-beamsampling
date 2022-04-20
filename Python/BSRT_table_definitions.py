### Beam Sampling RT tables

# RT model
import numpy as np
import random


from asrt_table_definitions import *
from rt_model import Seed, DataFilter
import beam_sampling_rt as bs_rt

schema = dj.schema('asrt', locals())

@schema
class BSRT_Hypers(dj.Lookup):
    definition = """
    hyper_id: int unsigned # value of the seed
    ---
    alpha0_a: double
    alpha0_b: double
    gamma_a: double
    gamma_b: double
    numb: int # number of burn-in steps
    nums: int # number of samples
    numi: int # number of steps between samples
    numi_rt_phi: int # number of steps for resampling rt and phi params
    tau0_a: double
    tau0_b: double
    mu_a: double
    mu_b: double
    sigma_a: double
    sigma_b: double
    h: blob  # prior for emissions
    k0: int unsigned
    """

@schema
class BSRT_Samples(dj.Computed):
    definition = """
    -> BSRT_Hypers
    -> Seed
    -> DataFilter.Dataset
    """

    class Sample(dj.Part):
        definition = """
        ->BSRT_Samples
        iteration: int unsigned # id of iteration
        ---
        loglike: double # log-likelihood
        tau0: double
        mu: double
        sigma: double
        beta: longblob
        alpha0: double
        gamma: double
        pi: longblob # transition matrix
        phi: longblob # emission matrix
        s: longblob # hidden state sequence
        """

    def _make_tuples(self, key):
        print('Populating  seed: '+str(key['seed'])+'   filter_id: '+str(key['filter_id'])+'   dataset_id: '+str(key['dataset_id']))
        self.insert1(key)
        seed = (Seed() & key).fetch1(as_dict=True)['seed']
        filter_details = (DataFilter() & key).fetch1(as_dict = True)
        hypers = (BSRT_Hypers() & key).fetch1(as_dict = True)
        hypers['H'] = hypers['h']
        hypers['rt_max'] = filter_details['rt_max']

        np.random.seed(seed)


        for data in DataFilter.Dataset.GetData(key):
            if 'NAPEXP' in data['participant_id'].iloc[0]:
                return()
            rt = np.array(data['first_rt'])
            Y  = np.array(data['event'])-1
            K  = hypers['k0']

            max_tries = 40
            successful = False
            tries = 0
            while tries < max_tries and not successful:
                try: 
                    S0 = np.array([np.random.choice(range(K)) for i in range(len(Y))])

                    samples, stats = bs_rt.sample_beam_ihmm(S0.copy(), Y, rt, hypers, hypers['numb'], hypers['nums'], hypers['numi'], hypers['numi_rt_phi'], hypers['numi_rt_phi'])
                    successful = True
                except Exception as e:
                    print('Unsuccessful sampling.')
                    print(e)
                    tries += 1

        if successful:
            entries = []
            for s in samples:
                entry = dict(beta = s['Beta'],
                             phi = s['Phi'],
                             pi = s['Pi'],
                             tau0 = s['tau0'],
                             mu = s['mu'],
                             sigma = s['sigma'],
                             iteration = s['iteration'],
                             loglike = s['loglike'],
                             s = s['S'],
                             gamma = s['gamma'],
                             alpha0 = s['alpha0'])
                entry = {**entry, **key}
                entries.append(entry)

            self.Sample().insert(entries)
