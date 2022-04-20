import numpy as np
import pickle
from definitions import OUTPUT_DIR
from os.path import join as pj
import glob

keys = ['pi', 'phi',
        'tau0', 'mu','sigma',
        'jll', 'K',
       ]

def faulty_chain(samples):
    if np.max(samples['tau0']) > 100:
        return(True)
    if np.max(samples['mu']) > 100:
        return(True)
    if np.max(samples['sigma']) > 100:
        return(True)
    if np.min(samples['mu']) < 0.00001:
        return(True)
    if np.min(samples['tau0']) < 0.0001:
        return(True)
    return(False)

def sample_filename(ini, participant, epoch, commit = None):
    _ini = ini.split('.')[0] if '.ini' in ini else ini
    if '_markov' in _ini:
        return(pj(commit,'_'.join([commit, _ini.split('_')[0], participant, 'blocks', epoch, 'markov_samples.pkl'])))
    if '_hmm' in _ini:
        return(pj(commit,'_'.join([commit, _ini.split('_')[0], participant, 'blocks', epoch, 'hmm'+_ini.split('hmm')[-1], 'samples.pkl'])))
    else:
        return(pj(commit,'_'.join([commit, _ini,participant,'blocks',epoch,'samples.pkl'])))

def sample_filename_list(inis=[], commits=[], participants=[], epochs=[]):
    l = glob.glob(OUTPUT_DIR+"/*")
    criteria = [inis,commits,participants,epochs]
    for criterium_sets in criteria:
        l = [item for item in l if any(criterium in item for criterium in criterium_sets)]
    return l

def import_samples(**kwargs):
    res = []

    keys = ['pi', 'phi',
        'tau0', 'mu','sigma',
        'jll', 'K',
       ]

    if 'filenames' not in kwargs:
        if ('ini' in kwargs) and ('participant' in kwargs) and ('epoch' in kwargs):
            if 'commit' in kwargs:
                samples_filenames = [sample_filename(kwargs['ini'], kwargs['participant'], kwargs['epoch'], kwargs['commit'])]
            else:
                samples_filenames = [sample_filename(kwargs['ini'], kwargs['participant'], kwargs['epoch'])]
        elif ('inis' in kwargs) and ('participants' in kwargs) and ('epochs' in kwargs):
            samples_filenames = [sample_filename(i, p, b) for i,p,b in zip(kwargs['inis'], kwargs['participants'], kwargs['epochs'])]
        else:
            raise NotImplemented('Use samples_filenames or ini, participant and epoch, or inis, participants, epochs arguments.')
    else:
        samples_filenames = kwargs['filenames']

    if 'keys' in kwargs:
        keys = kwargs['keys']

    for samples_filename in samples_filenames:
        with open(pj(OUTPUT_DIR,samples_filename),'rb') as file:
            res.append(pickle.load(file))
            samples = dict()
            chains = 0
            for key in keys:
                if np.all(np.array([(key in r['samples'][0]) for r in res])):
                    samples[key] = []

    for i in range(len(res)):
        _samples = dict()
        for key in samples.keys():
            _samples[key] = [s[key] for s in res[i]['samples']]
        if True: #not faulty_chain(_samples) and len(res[i]['samples'])==len(res[i]['stats']['jll']):
            for key in samples.keys():
                samples[key].extend(_samples[key])
            chains += 1
        else:
            print('Faulty chain: ', samples_filenames[i])

    if ('put_together' in kwargs) and (kwargs['put_together'] and (len(res)>1)):
        for r in res[1:]:
            res[0]['samples'].extend(r['samples'])
            res[0]['stats'].extend(r['stats'])

    if 'top_n_samples' in kwargs:
        top_n_samples = kwargs['top_n_samples']
        samples = np.array([])
        for c in range(chains):
            _samples = unique_samples(np.array(res[c]['samples']))
            _samples = _samples[np.argsort(res[c]['stats']['jll'])[-top_n_samples:]]
            samples = np.concatenate([samples, _samples])

    elif 'last_n_samples' in kwargs:
        samples = np.array([])
        for c in range(chains):
            _samples = unique_samples(np.array(res[c]['samples']))
            _samples = _samples[-kwargs['last_n_samples']:]
            samples = np.concatenate([samples, _samples])

    return(samples)

def top_jll_samples(filenames, n_samples=None):
    res, _, chains = import_samples(filenames=filenames)
    samples = np.array([])
    for c in range(chains):
        _samples = np.array(res[c]['samples'])
        if n_samples is not None:
            _samples = _samples[np.argsort(res[c]['stats']['jll'])[-n_samples:]]
        samples = np.concatenate([samples, _samples])
    return(samples)

def unique_samples(samples):
    tau0 = np.array([sample['tau0'] for sample in samples])
    _, indices = np.unique(tau0, return_index=True)
    filtered_samples = np.array(samples)[sorted(indices)]
    return(filtered_samples)
