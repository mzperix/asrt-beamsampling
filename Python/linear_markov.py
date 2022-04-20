import pickle
from definitions import OUTPUT_DIR
from os.path import join as pj
import numpy as np

class LinearMarkovSamples:
    def __init__(self, **kwargs):
        if 'samples_filenames' in kwargs:
            samples_filenames = kwargs['samples_filenames']
            self.load_from_samples_filenames(samples_filenames)
        
        if (('participant' in kwargs) 
            and ('ini' in kwargs) 
            and ('commit' in kwargs)
            and ('epoch' in kwargs)):
            
            participant = kwargs['participant']
            ini = kwargs['ini']
            commit = kwargs['commit']
            epoch = kwargs['epoch']
            samples_filename = commit+'_'+ini.split('.')[0]+'_'+participant+'_blocks_'+epoch+'_linear_markov_samples.pkl'
            self.load_from_samples_filenames([samples_filename])
        
        if 'last_n_samples' in kwargs:
            self.samples = self.samples[-kwargs['last_n_samples']:]
            
    def load_from_samples_filenames(self, samples_filenames):
        self.samples=[]
        for filename in samples_filenames:
            with open(pj(OUTPUT_DIR,filename),'rb') as file:
                self.samples.extend(pickle.load(file)['samples'])
        
    def instantiate_sample(self, sample):
        self.model = LinearMarkov(sample['mu'], sample['sigma'], sample['pi'])
        
    def predict_rt(self, Y):
        predictions = []
        for sample in self.samples:
            self.instantiate_sample(sample)
            predictions.append(self.model.predict_rt(Y))
        return(np.mean(np.array(predictions),0))
            
class LinearMarkov:
    def __init__(self, mu, sigma, pi):
        self.mu = mu
        self.sigma = sigma
        self.pi = pi
        
    def predict_rt(self, Y):
        rt = [self.mu+np.mean(self.pi[:,Y[0]])]
        y_previous = Y[0]
        for y in Y[1:]:
            rt.append(self.mu+self.pi[y_previous,y])
            y_previous = y
        return(np.array(rt))