import ideal_observer as i_o
import linear_markov
import triplets
import data_handling as dh
import pandas as pd
import numpy as np
from utils import pd_multiplicator, get_git_hash
import filters
import traceback
from definitions import OUTPUT_DIR
import argparse
from os.path import join as pj
from functools import partial, reduce
import pprint
from ideal_observer import STATE_PRIOR_ENTROPIES, STATE_POSTERIOR_ENTROPIES,\
    STATE_POSTERIOR_DIFFERENCE, PRED_DIFFERENCE, PRED_ENTROPIES, \
    PRED_PROB_ALL, PRED_PROBS, LOG_PRED_PROBS, LOG_PRED_PROB
import itertools
import scipy.stats as sp

MODELS = {'iHMM': i_o.IdealObserverSamples,
          'Markov': i_o.IdealObserverSamples,
          'LinearMarkov': linear_markov.LinearMarkovSamples,
          'Triplet': triplets.TripletModel,
          'GroundTruth': partial(i_o.GroundTruth,
              transition_noise=0.13,
              emission_noise=0.1)}

def data_and_model(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier, permute):

    if participant_test is None:
        participant_test = participant_train
    data_train = dh.import_data(participant = participant_train, block = e_train)
    data_test  = dh.import_data(participant = participant_test, block = e_test)

    if participant_test and permute:
        data_test = permute_keys(data_train, data_test)
    else:
        data_test['permuted'] = False

    if '.ini' in ini:
        ini = ini.split('.')[0]
    data_test['ini'] = ini
    data_test['model'] = model_name
    data_test['filters'] = filters.get_filter(data_test,ini)
    data_train['filters'] = filters.get_filter(data_train,ini)
    data_test['commit'] = commit
    data_test['e_train'] = e_train
    data_test['e_test'] = e_test
    data_test['participant_test'] = participant_test.split('_')[-1]
    data_test['participant_train'] = participant_train.split('_')[-1]

    if model_name == 'Markov':
        ini = ini+'_markov'

    kwargs = {}
    if model_name == 'GroundTruth':
        kwargs['sequence'] = data_train['Y'][data_train['trial_type']=='P'][:4]
    model = MODELS[model_name](participant = modifier+participant_train,
        ini = ini, commit = commit, epoch = e_train, last_n_samples=last_n_samples,
        **kwargs)
    return data_train, data_test, model

@pd_multiplicator
def residuals(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier='', permute=False):
    data_train, data_test, model = data_and_model(
            participant_train, participant_test,
            ini, commit, e_train, e_test, model_name,
            last_n_samples, modifier, permute)
    data_test['r'] = model.calculate_residuals(data_test)
    data_test['r'] = np.mean(data_test['r'],0)
    df = pd.DataFrame(data_test)
    return df

@pd_multiplicator
def entropy(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier='', permute=False):
    assert (model_name == 'iHMM' or model_name=='GroundTruth')
    try:
        data_train, data_test, model = data_and_model(
            participant_train, participant_test,
            ini, commit, e_train, e_test, model_name,
            last_n_samples, modifier, permute)

        prediction = i_o.state_and_prediction_entropy(model, data_test['Y'])
        for key in [STATE_POSTERIOR_ENTROPIES, STATE_PRIOR_ENTROPIES,\
            STATE_POSTERIOR_DIFFERENCE, PRED_ENTROPIES, PRED_DIFFERENCE,\
            PRED_PROBS]:
            data_test[key] = np.mean(prediction[key],0)

        data_test[PRED_ENTROPIES] = prediction[PRED_ENTROPIES]
        data_test['K_mean'] = np.mean([np.size(sample['phi'],0) for sample in model.samples])
        df = pd.DataFrame(data_test)

    except Exception as e:
        print('Failed: p_train', participant_train,'p_test', participant_test,
            'i', ini, 'c', commit, 'e_train', e_train, 'e_test', e_test,
            'model_name', model_name, 'last_n_samples', last_n_samples)
        print(e)
        tb = traceback.format_exc()
        print(tb)
        df = pd.DataFrame()
    return df

@pd_multiplicator
def predicted_rt(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier='', permute=False):

    try:
        data_train, data_test, model = data_and_model(
            participant_train, participant_test,
            ini, commit, e_train, e_test, model_name,
            last_n_samples, modifier, permute)

        if model_name == 'Triplet':
            data_test['rt_predicted'] = model.predict_rt(data_test, data_train)
        else:
            data_test['rt_predicted'] = model.predict_rt(data_test['Y'])
        df = pd.DataFrame(data_test)

    except Exception as e:
        print('Failed: p_train', participant_train,'p_test', participant_test,
            'i', ini, 'c', commit, 'e_train', e_train, 'e_test', e_test,
            'model_name', model_name, 'last_n_samples', last_n_samples)
        print(e)
        tb = traceback.format_exc()
        df = pd.DataFrame()
        print(tb)
    return df

@pd_multiplicator
def fingerprint(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier='', permute=False, length=3):

    def add_sequence_to_data(data, length):
        seq = ''
        data['sequence'] = ['']*len(data['Y'])
        for i,y in enumerate(data['Y']):
            # RESTART STREAM MISSING BUT SEQ LENGTHS ARE SHORTER THAN FILTERED
            # STARTS OF EACH BLOCK
            if len(seq)==length:
                seq = seq[1:]
            seq += str(y+1)
            data['sequence'][i] = seq
        return data

    try:
        data_train, data_test, model = data_and_model(
            participant_train, participant_test,
            ini, commit, e_train, e_test, model_name,
            last_n_samples, modifier, permute)

        sequences = [seq for seq in itertools.product(*[[0,1,2,3]]*length)]
        data = {key: data_test[key] \
            for key in ['commit', 'participant_test', 'model', 'e_train', 'ini']
        }

        data_test = add_sequence_to_data(data_test, length)
        df = pd.DataFrame(data_test)
        df = df[df.correct_response==1]
        df = df[df.filters==True]
        df = df.groupby(['sequence'])\
            [["rt"]].agg([np.mean,sp.sem,'count']).reset_index()
        df.columns = ["_".join(x) for x in df.columns.ravel()]

        #print(df.head())
        data['sequence'] = []
        data['predicted_probability'] = []
        data['log_predicted_probability'] = []
        data['rt_predicted'] = []

        for seq in sequences:
            data['sequence'].append(
                ''.join(map(lambda x: str(x+1),seq)))
            prediction = model.seed_and_predict(seq)
            data['predicted_probability'].append(np.mean(prediction['predicted_probability']))
            data['log_predicted_probability'].append(np.mean(np.log(prediction['predicted_probability'])))
            data['rt_predicted'].append(np.mean(prediction['rt']))
        df = pd.merge(pd.DataFrame(data),df,
            left_on='sequence', right_on='sequence_', how='outer')

    except Exception as e:
        print('Failed: p_train', participant_train,'p_test', participant_test,
            'i', ini, 'c', commit, 'e_train', e_train, 'e_test', e_test,
            'model_name', model_name, 'last_n_samples', last_n_samples)
        print(e)
        tb = traceback.format_exc()
        df = pd.DataFrame()
        print(tb)
    return df

def rt_correlation(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier='', permute=False, parallel=1):

    df = predicted_rt(participant_train, participant_test,
        ini, commit, e_train, e_test, model_name,
        last_n_samples, modifier, permute, parallel=parallel)
    df = df[df.filters==True]
    df = df[df.correct_response==1]
    print('Computing correlations.')
    df = df.groupby(['participant_train','participant_test',
        'commit','e_train','e_test','model','permuted','block','trial'])\
        [["rt", "rt_predicted"]].mean().reset_index()

    # CUTOFF LARGE RTs
    df['z'] = df.groupby(
        ['participant_train','participant_test',
        'commit','e_train','e_test','model','permuted'],
        group_keys=False)\
        .apply(lambda g: (g.rt-g.rt.mean())/g.rt.std())
    df = df[df.z < 3.0]

    del df['z']

    df = df.groupby(['participant_train','participant_test',
        'commit','e_train','e_test','model','permuted'])\
        [["rt","rt_predicted"]]\
        .corr('pearson')[["rt_predicted"]].iloc[0::2,:]\
        .rename(columns={"rt_predicted": "correlation"})
    print(df.head())
    return df

@pd_multiplicator
def pred_probs(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier='', permute=False):

    try:
        data_train, data_test, model = data_and_model(
            participant_train, participant_test,
            ini, commit, e_train, e_test, model_name,
            last_n_samples, modifier, permute)
        if model_name == 'Triplet':
            predictions = model.generate_predictive_probabilities(data_test)
        else:
            predictions = model.generate_predictive_probabilities(data_test['Y'])
            data_test[LOG_PRED_PROB] = np.mean(predictions[LOG_PRED_PROBS],0)
        pred_prob_all = np.array(predictions[PRED_PROB_ALL])
        for event in range(4):
            data_test['y'+str(event)] = np.mean(pred_prob_all,0)[event,:]
        df = pd.DataFrame(data_test)

    except Exception as e:
        print('Failed: p_train', participant_train,'p_test', participant_test,
            'i', ini, 'c', commit, 'e_train', e_train, 'e_test', e_test,
            'model_name', model_name, 'last_n_samples', last_n_samples)
        print(e)
        tb = traceback.format_exc()
        df = pd.DataFrame()
        #print(tb)
    return df

@pd_multiplicator
def model_samples(participant_train, participant_test,
    ini, commit, e_train, e_test, model_name,
    last_n_samples, modifier='', permute=False, n_models=10):

    try:
        data_train, data_test, model = data_and_model(
            participant_train, participant_test,
            ini, commit, e_train, e_test, model_name,
            last_n_samples, modifier, permute)

        d = {key: [sample[key] for sample in model.samples[-n_models:]]
            for key in model.samples[-1].keys()}
        if not ('iteration' in d):
            d['iteration'] = range(len(model.samples)-n_models,len(model.samples))
        d['participant_train'] = participant_train
        d['ini'] = ini
        d['e_train'] = e_train
        d['model'] = model_name
        d['commit'] = commit
        d['sequence'] = [(data_train['Y'][data_train['trial_type']=='P'])[0:4]] * n_models
        #print(d)
        df = pd.DataFrame(d)

    except Exception as e:
        print('Failed: p_train', participant_train,'p_test', participant_test,
            'i', ini, 'c', commit, 'e_train', e_train, 'e_test', e_test,
            'model_name', model_name, 'last_n_samples', last_n_samples)
        #print(e)
        tb = traceback.format_exc()
        df = pd.DataFrame()
        print(tb)
    return df

def permute_keys(data_template,data):
    def keymap(key,sequence_template, sequence):
        return sequence_template[np.where(sequence==key)][0]

    sequence_template = np.array(data_template['Y'][data_template['trial_type']=='P'][:5])
    sequence = np.array(data['Y'][data['trial_type']=='P'][:5])
    _data = data.copy()
    _data['Y'] = np.array([keymap(k, sequence_template, sequence) for k in data['Y']])
    _data['permuted'] = True
    return _data

def prediction_correlation(df):
    return(df[(df.filters==True) & (df.correct_response==1)].groupby( \
           ['participant','ini','commit','e_train','model'])[['rt','rt_predicted']].corr( \
           )['rt_predicted'][0::2].reset_index().rename(columns={'rt_predicted': 'performance'}))

def map_response(r):
    if (r == 1) or (r == 'z'):
        return(0)
    if (r == 2) or (r == 'c'):
        return(1)
    if (r == 4) or (r == 'b'):
        return(2)
    if (r == 5) or (r == 'm'):
        return(3)

def error_performance(df):
    prediction_histogram = []
    predictions = np.vstack([df['y'+str(i)] for i in range(4)])
    for t in range(len(df)):
        Y = df.Y.iloc[t]
        prediction_histogram.append(np.where(np.sort(predictions[:,t])[::-1]==predictions[map_response(df['first_response'].iloc[t]),t])[0][0]+1)
    df['model_event'] = prediction_histogram
    return(df)

def compute_results(commit, template_name, parallel):
    """Generates all results csv files for the figures and statistics.
       Inputs:
            commit: first 7 character id of commit
            template_name: result computation parameter settings
            parallel: number of parallel threads to use"""

    FUNCTIONS = {'predicted_rt': predicted_rt,
                 'predicted_prob': pred_probs,
                 'residuals': residuals,
                 'entropy': entropy,
                 'rt_correlation': rt_correlation,
                 'fingerprint': fingerprint,
                 'model_samples': partial(model_samples, n_models=5)}

    template = dh.get_template(template_name)
    experiment = template['EXPERIMENT']['name']

    if 'participant_test' in template['EXPERIMENT']:
        if template['EXPERIMENT']['participant_test'].isdigit():
            participant_test = ['_'.join([experiment,p]) for p in template['EXPERIMENT']['participant_test'].split(',')]
        else:
            participant_test = template['EXPERIMENT']['participant_test'].split(',')
    else:
        participant_test = None

    if 'participant_train' in template['EXPERIMENT']:
        participant_train = ['_'.join([experiment,p]) for p in template['EXPERIMENT']['participant_train'].split(',')]
    else:
        participant_train = ['_'.join([experiment,p]) for p in template['EXPERIMENT']['participants'].split(',')]

    last_n_samples = int(template['SAMPLES']['last_n_samples'])
    inis = template['SAMPLES']['inis'].split(',')
    permute = False

    git_hash = get_git_hash()
    result_names = template.keys()
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(template)

    for result_name in result_names:
        if not result_name in ['EXPERIMENT','SAMPLES','DEFAULT']:
            # PARTICIPANTS
            #participant_train = participants
            if template[result_name]['design_subject'] == 'across':
                if participant_test is None:
                    participant_test = participant_train
                permute=[]
                if 'T' in template[result_name]['permute']:
                    permute.append(True)
                if 'F' in template[result_name]['permute']:
                    permute.append(False)
            else:
                participant_test = None

            # FUNCTION
            fun = partial(FUNCTIONS[template[result_name]['function']], parallel=parallel)
            if template[result_name]['function'] == 'fingerprint':
                fun = partial(fun,length=int(template[result_name]['length']))
            model_names  = template[result_name]['models'].split(',')
            epoch_trains = template[result_name]['epoch_train'].split(',')
            epoch_tests  = template[result_name]['epoch_test'].split(',')

            if template[result_name]['design_session'] == 'within':
                df = pd.DataFrame()
                for e_train, e_test in zip(epoch_trains, epoch_tests):
                    df = pd.concat([df,fun(participant_train, participant_test, inis, commit,
                        e_train, e_test, model_names, last_n_samples, permute=permute)])
            if template[result_name]['design_session'] == 'across':
                df = fun(participant_train, participant_test, inis, commit,
                    epoch_trains, epoch_tests, model_names, last_n_samples, permute=permute)

            print('Finished', result_name)
            output = pj(OUTPUT_DIR,commit,'_'.join([git_hash, result_name, template_name]))
            if ('file_type' in template[result_name]):
                file_type = template[result_name]['file_type']
                output += '.' + file_type
                if file_type == 'json':
                    df.to_json(output, orient='records')
                    print('Saved results in', output)
                else:
                    print('UNKNOWN FILE FORMAT.')
            else:
                output += '.csv'
                df.to_csv(output, float_format="%g")
                print('Saved results in', output)

def main(args):
    compute_results(args.c, args.t, args.p)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Compute data for plots and statistics by template.')
    parser.add_argument('-c', type=str, default=None, help="Commit id (first 7 characters)")
    parser.add_argument('-t', type=str, default=None, help="Template filename")
    parser.add_argument('-p', type=int, default=1, help="Number of parallel processes")
    args = parser.parse_args()
    main(args)
