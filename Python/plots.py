import sys
sys.path.append('../../')
import stan_plots as s_p
from matplotlib import rcParams
import ideal_observer as i_o
import pickle
import numpy as np
import matplotlib.pyplot as plt
from triplets import triplet_model_predict, triplet_correlation
import scipy.stats as sp

colors = ['#34a5da', '#84b7c0', '#fdd262', '#df422a']
SMALL_SIZE = 10
MEDIUM_SIZE = 14
BIGGER_SIZE = 20

N_MODEL_SAMPLES = 40

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def QQ_plot(samples1, samples2, quantiles=[0.05,0.15,0.25,0.35,0.42,0.5,0.58,0.65,0.75,0.85, 0.95], show = True):
	x, y = [], []
	_s1 = np.sort(samples1)
	_s2 = np.sort(samples2)
	l1 = len(_s1)
	l2 = len(_s2)
	for q in quantiles:
		x.append(_s1[int(np.floor(l1*q))])
		y.append(_s2[int(np.floor(l2*q))])
	x0 = np.min(np.append(x,y))
	x1 = np.max(np.append(x,y))
	plt.plot([x0,x1],[x0,x1])
	plt.plot(x,y)
	plt.plot(x,y, '.')
	if(show):
		plt.show()
        
def norm_rows(m):
    return((m.transpose()/np.sum(m,1)).transpose())

def cutoff_matrix(m, epsilon):
    return(norm_rows(m*np.array(m>epsilon)))

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

def stacked_bar_plot(ax, pred_probs):
    T = np.size(pred_probs,1)
    x = np.arange(T)
    current_means = np.repeat(0.0,T)
    for i in range(4):
        ax.bar(x, pred_probs[i,:],
               bottom = current_means,
               width = 1,
               color = colors[i])
        current_means += pred_probs[i,:]

def dot_plot(ax, pred_probs, data, direction = 'horizontal'):
    T = np.size(pred_probs,1)
    Y = data['Y'][:T]
    pred_probs -= np.min(pred_probs)
    pred_probs /= np.max(pred_probs)
    #ax.patch.set_facecolor('#000000')
    if direction == 'horizontal':
        x = np.arange(T)
        for i in range(np.max(Y)+1):
            ax.scatter(x, 
                       np.repeat(i,T),
                      # s = 200+30*np.log(pred_probs[i,:]), 
                       s = 150*pred_probs[i,:], 
                       c = colors[i])
        values = np.array([pred_probs[Y[t],t] for t in np.where(data['trial_type'][:T] == 'P')])

        ax.scatter(x[data['trial_type'][:T] == 'P'], 
                   Y[data['trial_type'][:T] == 'P'],
                   s = values*52, 
                   c = '#ffffff')
        values = np.array([pred_probs[Y[t],t] for t in np.where(data['trial_type'][:T] == 'R')])
        ax.scatter(x[data['trial_type'][:T] == 'R'], 
                   Y[data['trial_type'][:T] == 'R'],
                   s = values*52,
                   marker = 'x',
                   c = '#ffffff')
        ax.set_ylim(-0.5,3.5)
        
    if direction == 'vertical':
        y = np.arange(T)[::-1]
        for i in range(np.max(Y)+1):
            ax.scatter(np.repeat(i,T), 
                       y,
                       s = 150*pred_probs[i,:], 
                       c = colors[i])
        values = np.array([pred_probs[Y[t],t] for t in np.where(data['trial_type'][:T] == 'R')[0]])
        ax.scatter(Y[data['trial_type'][:T] == 'R'],
                   y[data['trial_type'][:T] == 'R'],
                   s = values*52, 
                   c = '#ffffff')
        ax.set_xlim(-0.5,3.5)
    
keys = ['pi', 'phi', 
        'tau0', 'mu','sigma',
        'jll', 'K', 
       ]

def import_samples(samples_filenames, keys = keys, put_together=False):
    res = []
    for samples_filename in samples_filenames:
        with open(samples_filename,'rb') as file:
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
        if not faulty_chain(_samples) and len(res[i]['samples'])==len(res[i]['stats']['jll']):
            for key in samples.keys():
                samples[key].extend(_samples[key])
            chains += 1
        else:
            print('Faulty chain: ', samples_filenames[i])
    if put_together and len(res)>1:
        for r in res[1:]:
            res[0]['samples'].extend(r['samples'])
            res[0]['stats'].extend(r['stats'])
    return(res, samples, chains)

def stan_plots(samples, chains, keys):
    ## STAN RESULT PLOTS
    stan_plot_samples = {k: samples[k] for k in samples.keys() & keys}
    s_p.sample_plots(samples = stan_plot_samples, 
                     chains = chains,  
                     plots = ['bar','chains','scatter_matrix'])
    
def model_plots(res, chains, data, test_data = None, prediction = 'MARGINAL'):
    incorrect_means = []
    incorrect_stds = []
    correct_means = []
    correct_stds = []

    test_correlations, triplet_correlations = [], []
    
    width = 1 # Train correlations
    if prediction == 'MAP':
        width += 2 # Matrix plots
    if test_data is not None:
        width += 2 # Test correlations plus miss predictions
    
    rcParams['figure.figsize'] = 5*width, chains*5
    fig, ax = plt.subplots(chains, width)
    ax = ax.reshape((chains, width)) # if chains == 1
    fig.subplots_adjust(hspace=0.32, wspace=0.22)

    Y = data['Y']
    T = np.size(Y)
    rt = data['rt'][data['correct_response']==1]
    if test_data is not None:
        test_Y = test_data['Y']
        test_T = np.size(Y)
        test_rt = test_data['rt'][test_data['correct_response']==1]
    
    for c in range(chains):
        if prediction == 'MAP':
            samples = res[c]['samples'][np.argmax(res[c]['stats']['jll'])]
            pi_mean = samples['Pi']
            ax[c,-2].matshow(pi_mean,
                            vmin = 0, vmax = 1)
            ax[c,-2].set_title('pi')
            phi_mean = samples['Phi']
            ax[c,-1].matshow(phi_mean,
                            vmin = 0, vmax = 1)
            ax[c,-1].set_title('phi')
            samples = [samples]
            
        if prediction == 'MARGINAL':
            samples = np.array(res[c]['samples'])
            samples = samples[np.argsort(res[c]['stats']['jll'])[-N_MODEL_SAMPLES:]]
        
        sample_log_pred_prob = i_o.marginal_prediction(samples, Y)
        sample_log_pred_prob = sample_log_pred_prob[data['correct_response']==1]        
        ax[c,0].text(0.06,
                     0.94,
                     'Corr: %.4f' % np.corrcoef(sample_log_pred_prob, rt)[0,1],
                     transform = ax[c,0].transAxes,
                     size = MEDIUM_SIZE)

        ax[c,0].plot(sample_log_pred_prob, rt,'.', markeredgecolor = colors[0], markerfacecolor = colors[0])
        ax[c,0].set_title('TRAIN')
        ax[c,0].set_xlabel('predicted log probability')
        ax[c,0].set_ylabel('reaction time')
        ax[c,0].set_ylim(np.mean(rt)-3*np.std(rt), np.mean(rt)+3*np.std(rt))
        if test_data is not None:
            test_Y = test_data['Y']
            test_sample_log_pred_prob = i_o.marginal_prediction(samples, test_Y)
            test_sample_log_pred_prob = test_sample_log_pred_prob[test_data['correct_response']==1]
            ax[c,1].text(0.06,
                         0.94,
                         'Corr: %.4f' % np.corrcoef(test_sample_log_pred_prob, test_rt)[0,1],
                         transform = ax[c,1].transAxes,
                         size = MEDIUM_SIZE)
            test_correlations.append(np.corrcoef(test_sample_log_pred_prob, test_rt)[0,1])
            ax[c,1].text(0.06,
                         0.88,
                         'Triplet corr: %.4f' % triplet_correlation(test_data),
                         transform = ax[c,1].transAxes,
                         size = MEDIUM_SIZE)
            triplet_correlations.append(triplet_correlation(test_data))
            ax[c,1].plot(test_sample_log_pred_prob, test_rt,'.', markeredgecolor = colors[0], markerfacecolor = colors[0])
            ax[c,1].set_title('TEST')
            ax[c,1].set_xlabel('predicted log probability')
            ax[c,1].set_ylabel('reaction time')
            ax[c,1].set_ylim(np.mean(test_rt)-3*np.std(test_rt), np.mean(test_rt)+3*np.std(test_rt))

            ## Correct vs. incorrect predicted probabilities
            test_sample_log_pred_prob = i_o.marginal_prediction(samples, test_Y)
            pred_incorrect = test_sample_log_pred_prob[(test_data['correct_response']==0)*(test_data['trial_type']=='R')]
            pred_correct = test_sample_log_pred_prob[(test_data['correct_response']==1)*(test_data['trial_type']=='R')]
            ax[c,2].boxplot([pred_incorrect, pred_correct], showmeans = True)
            ax[c,2].set_title('Incorrect vs. correct pred prob')
            ax[c,2].text(0.06,
                         0.94,
                         'No. incorrect trials: %i' % len(pred_incorrect),
                         transform = ax[c,2].transAxes,
                         size = MEDIUM_SIZE)
            
            incorrect_means.append(np.mean(pred_incorrect))
            incorrect_stds.append(np.std(pred_incorrect)/np.sqrt(np.sum(test_data['correct_response']==0)))
            correct_means.append(np.mean(pred_correct))
            correct_stds.append(np.std(pred_correct)/np.sqrt(np.sum(test_data['correct_response']==1)))
    plt.show()
    return((incorrect_means, incorrect_stds, correct_means, correct_stds,
            test_correlations, triplet_correlations))

def raw_data_plot(data, ax = None, show = True, title = ''):
    if ax is None:
        fig, ax = plt.subplots()
        rcParams['figure.figsize'] = 8,6
    T = len(data['rt'])
    T = 200
    ax.plot(np.where(data['correct_response'][:T]==1)[0], 
            data['rt'][:T][data['correct_response'][:T]==1], 
            '.',
            markeredgecolor = colors[0], 
            markerfacecolor = colors[0],
            markersize = 8,
           )
    b = min(data['rt'][:T])
    #ax.scatter(x = range(T), 
    #           y = [b-40+y*10 for y in data['Y'][:T]], 
    #           s = 20,
    #           marker = 's',
    #           color = [colors[y] for y in data['Y'][:T]],
    #          )
    #ax.scatter(x = range(T), 
    #           y = data['rt'][:T], 
    #           s = 12,
    #           color = [colors[y] for y in data['Y'][:T]],
    #          )
    
    
    
    ax.set_title(title)
    ax.set_xlabel('Trial')
    ax.set_ylabel('Reaction time')
    if show:
        plt.show()

def prediction_plot(res, chains, data, T, prediction = 'MARGINAL'):
    rcParams['figure.figsize'] = 20*T/70, 1.4+1.4*chains    
    fig, ax = plt.subplots(chains+1,1)
    Y = data['Y']
    stimuli = np.zeros((np.max(Y)+1,T))
    for t in range(T):
        stimuli[Y[t],t] = 1
    dot_plot(ax[0], stimuli, data)
    ax[0].axis('off')
    ax[0].set_title('Actual stimuli')
    
    for c in range(chains):
        if prediction == 'MAP':
            samples = [res[c]['samples'][np.argmax(res[c]['stats']['jll'])]]
        if prediction == 'MARGINAL':
            samples = np.array(res[c]['samples'])
            samples = samples[np.argsort(res[c]['stats']['jll'])[-N_MODEL_SAMPLES:]]
        
        log_pred_probs = i_o.marginal_prediction_full(samples, Y[:T])
        dot_plot(ax[c+1], np.exp(log_pred_probs), data)
        #fig.patch.set_visible(False)
        ax[c+1].axis('off')
    plt.show()
        
def plot_train_result(data_filename, samples_filenames, test_data_filename = None, keys = keys, diagnostics = True, prediction = 'MARGINAL', n_show_chains = 0):    
    # LOAD DATA
    with open(data_filename,'rb') as file:
        data = pickle.load(file)
    
    if test_data_filename is not None:
        with open(test_data_filename, 'rb') as file:
            test_data = pickle.load(file)
    
    res, samples, chains = import_samples(samples_filenames, keys)
    print('Chains:', chains)
    print('len(res):', len(res))
    
    if test_data_filename is None:
        test_data = None
    
    if chains < 1:
        return()
    if n_show_chains > 0:
        chains = n_show_chains
    
    if diagnostics:
        stan_plots(samples, chains, keys = ['K','tau0','mu','sigma','jll','epsilon'])
        for c in range(chains):
            print('TOTAL SAMPLING TIME (chain',c,'): ', np.sum(res[c]['stats']['stan_time']), '  AVG. PER SAMPLE: ', np.sum(res[c]['stats']['stan_time'])/len(res[c]['samples']))
    
    # RAW DATA
    if test_data_filename is None:
        rcParams['figure.figsize'] = 6,4
        plot_raw_data(data, show = True, title = 'TRAIN')
    else:
        rcParams['figure.figsize'] = 12,4
        fig, axes = plt.subplots(1,2)
        raw_data_plot(data, axes[0], show = False, title = 'TRAIN')
        raw_data_plot(test_data, axes[1], show = True, title = 'TEST')
    
    print('NUMBER OF RT SAMPLES: ', len(data['rt']))
        
    print('Chain means')
    
    incorrect_means, incorrect_stds, correct_means, correct_stds, test_correlations, triplet_correlations = model_plots(res, chains, data, test_data, prediction)
    
    # PREDICTIONS WITH THE DOTPLOT
    prediction_plot(res, chains, test_data, T = 60)
    return((incorrect_means, incorrect_stds, correct_means, correct_stds,
            test_correlations, triplet_correlations))

def model_prediction_by_correctness_plot(incorrect_means, incorrect_stds, correct_means, correct_stds, pmin=-2.75, pmax=-1.1):
    rcParams['figure.figsize'] = 6,6
    plt.plot([pmin, pmax],[pmin, pmax])
    plt.errorbar(incorrect_means, correct_means, 2*np.array(incorrect_stds), 2*np.array(correct_stds), 
                 capsize=0, ls='none', color='#cccccc', 
                 elinewidth=1)
    plt.plot(incorrect_means, correct_means,'.', 
             markersize = 8, markeredgecolor = colors[0], markerfacecolor = colors[0])
    plt.xlabel('Incorrect mean log pred')
    plt.ylabel('Correct mean log pred')
    plt.title('Random trials', size = BIGGER_SIZE)
    plt.xlim(pmin,pmax)
    plt.ylim(pmin,pmax)
    plt.show()
    
def model_correlations(data_filenames, samples_filenames):
    data_list, samples_list = [], []
    for f in data_filenames:
        with open(f,'rb') as file:
            data_list.append(pickle.load(file))
    res, _, _ = import_samples(samples_filenames, keys = keys)
    N = len(data_list)
    correlations = np.zeros((N,N))
    for i, r in enumerate(res):
        samples = np.array(r['samples'])
        samples = samples[np.argsort(r['stats']['jll'])[-N_MODEL_SAMPLES:]]
        for j, data in enumerate(data_list):
            predictions = i_o.marginal_prediction(samples, data['Y'])
            correlations[i,j] = np.corrcoef(predictions[data['correct_response']==1], data['rt'][data['correct_response']==1])[0,1]
    rcParams['figure.figsize'] = 6,6
    plt.matshow(-correlations, vmin = 0, vmax = 1)
    plt.xlabel('Data')
    plt.ylabel('Samples')
    plt.title('Model correlations')
    plt.show()
    return(correlations)

def map_response(r):
    if r == 1:
        return(0)
    if r == 2:
        return(1)
    if r == 4:
        return(2)
    if r == 5:
        return(3)
    
def incorrect_predictions(data_filenames, samples_filenames, show = True):
    data_list, samples_list = [], []
    for f in data_filenames:
        with open(f,'rb') as file:
            data_list.append(pickle.load(file))
            #print(data_list[-1])
    res, _, _ = import_samples(samples_filenames, keys = keys)
    prediction_histogram = []
    for i, r in enumerate(res):
        samples = np.array(r['samples'])
        samples = samples[np.argsort(r['stats']['jll'])[-N_MODEL_SAMPLES:]]
        for j, data in enumerate(data_list):
            predictions = i_o.marginal_prediction_full(samples, data['Y'])
            responses = np.array([])
            for t, y in enumerate(data['Y']):
                if data['correct_response'][t] == 0:
                    prediction_histogram.append(np.where(np.sort(predictions[:,t])[::-1]==predictions[map_response(data['first_response'][t]),t])[0][0]+1)
    rcParams['figure.figsize']=6,4
    if show:
        plt.hist(prediction_histogram)
        plt.show()
    return(prediction_histogram)

def incorrect_predictions(data_filenames, samples_filenames, show = True):
    data_list, samples_list = [], []
    for f in data_filenames:
        with open(f,'rb') as file:
            data_list.append(pickle.load(file))
            #print(data_list[-1])
    res, _, _ = import_samples(samples_filenames, keys = keys)
    prediction_histogram = []
    for i, r in enumerate(res):
        samples = np.array(r['samples'])
        samples = samples[np.argsort(r['stats']['jll'])[-N_MODEL_SAMPLES:]]
        for j, data in enumerate(data_list):
            predictions = i_o.marginal_prediction_full(samples, data['Y'])
            responses = np.array([])
            for t, y in enumerate(data['Y']):
                if data['correct_response'][t] == 0:
                    predictions[data['Y'][t],t] = -np.inf
                    prediction_histogram.append(np.where(np.sort(predictions[:,t])[::-1]==predictions[map_response(data['first_response'][t]),t])[0][0]+1)
    rcParams['figure.figsize']=6,4
    if show:
        plt.hist(prediction_histogram)
        plt.show()
    return(prediction_histogram)

def incorrect_probabilities(data_filenames, samples_filenames):
    data_list, samples_list = [], []
    for f in data_filenames:
        with open(f,'rb') as file:
            data_list.append(pickle.load(file))
            #print(data_list[-1])
    res, _, _ = import_samples(samples_filenames, keys = keys)
    prediction_histogram = []
    for i, r in enumerate(res):
        samples = np.array(r['samples'])
        samples = samples[np.argsort(r['stats']['jll'])[-N_MODEL_SAMPLES:]]
        for j, data in enumerate(data_list):
            predictions = i_o.marginal_prediction_full(samples, data['Y'])
            responses = np.array([])
            for t, y in enumerate(data['Y']):
                if data['correct_response'][t] == 0:
                    prediction_histogram.append(predictions[map_response(data['first_response'][t]),t])
    prediction_histogram = np.exp(prediction_histogram)
    rcParams['figure.figsize']=6,4
    plt.hist(prediction_histogram)
    plt.show()
    return(prediction_histogram)


def latent_entropy(data_filenames, samples_filenames):
    # Take mean over state posteriors to get 
    # a mean representation weight and look at entropy of that
    data_list, samples_list = [], []
    for f in data_filenames:
        with open(f,'rb') as file:
            data_list.append(pickle.load(file))
            #print(data_list[-1])
    res, _, _ = import_samples(samples_filenames, keys = keys)
    entropies = []
    for i, r in enumerate(res):
        samples = np.array(r['samples'])
        samples = samples[np.argsort(r['stats']['jll'])[-5:]]
        for j, data in enumerate(data_list):
            entropies.append(i_o.filtered_state_posterior_entropy(samples, data['Y']))
            #plt.title(str(i))
            #plt.show()
            
    rcParams['figure.figsize']=2*len(entropies),4
    plt.boxplot(entropies)
    plt.show()
    return(entropies)