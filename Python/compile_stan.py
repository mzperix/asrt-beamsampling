from pystan import StanModel
import pickle
import sys

if len(sys.argv) > 1:
    stan_models = sys.argv[1:]
else:
    stan_models = ['sm_hypers', 'sm_pi_phi_rtparams', 'sm_pi_rtparams', 'sm_linear_markov']

for model in stan_models:
    print('Compiling model',model+'.stan')
    sm = StanModel(model+'.stan')
    with open(model+'.pkl','wb') as file:
        pickle.dump(sm, file)
    print('Stan model',model,' successfully pickled.')