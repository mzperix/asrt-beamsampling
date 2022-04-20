import data_handling as dh
import pandas as pd
from definitions import OUTPUT_DIR, DATA_DIR
from os.path import join as pj
import argparse
import glob
from utils import ton
import pickle
import traceback
import numpy as np
import samples_utils

def fnum(s):
    if s:
        if len(s)<=7:
            return s
        else:
            return "{:.1E}".format(float(s))
    else:
        return None

def parameter_settings(output=None):
    symbols = {'alpha0': u'\u03b1',
               'gamma': u'\u03b3',
               'beta_softmax': u'\u03b2'+u'soft',
               'numb': 'b',
               'nums': 's',
               'numi_pi_phi_rtparams': 'i',
               'new_stream_at_lowest_trial_in_block': 'cut',
               'error_sigma': u'\u03c3'+u'err',
               'tau0_shape': u'k\u03c4',
               'tau0_scale': u'\u03d1\u03c4',
               'mu_shape': u'k\u03bc',
               'mu_scale': u'\u03d1\u03bc',
               'sigma_shape': u'k\u03c3',
               'sigma_scale': u'\u03d1\u03c3',
               'skip_start_at_block': 'skip',
               'numi_hypers': 'i_hyp',
              }
    hypers_columns = ['alpha0', 'gamma', 
                      'h', 
                      'tau0_scale', 'tau0_shape', 'mu_scale','mu_shape','sigma_scale','sigma_shape',
                      'error_sigma',
                      'beta_softmax',
                      'new_stream_at_lowest_trial_in_block',
                      'rt_max']
    sampling_columns = ['k0', 'seed',
                        'numb', 'nums', 'numi_pi_phi_rtparams', 'numi_hypers']
    filters_columns = ['rt_min', 'skip_start_at_block']
    parameters = {}
    table_data = [['name'] + hypers_columns + sampling_columns + filters_columns]
    for i,item in enumerate(table_data[0]):
        if item in symbols.keys():
            table_data[0][i] = symbols[item]
  
    inis = glob.glob(pj(DATA_DIR,'*.ini'))
    for ini in inis:
        parameters = dh.get_ini(ini)
        row = [[ini.split('/')[-1].split('.')[0]] + 
               [fnum(ton(ton(parameters,'HYPERS'),col)) for col in hypers_columns] +
               [ton(ton(parameters,'SAMPLING_PARAMS'),col) for col in sampling_columns] +
               [ton(ton(parameters,'FILTERS'),col) for col in filters_columns]]
        table_data += row
    df = pd.DataFrame(table_data[1:], columns=table_data[0])
    if output:
        df.to_csv(pj(OUTPUT_DIR,output), float_format="%g", index=False)
        print('Generated parameter settings csv: '+pj(OUTPUT_DIR,output))
    return df

def sample_stats(commit, ignore=None):
    """Output:
          participant
          ini
          epoch
          parameter (tau0, sigma, mu, k, jll)
          mean
          sd: standard deviation
          q05: 5th percentile
          q95: 95th percentile"""
    files = glob.glob(pj(OUTPUT_DIR,commit,'*samples.pkl'))
    if ignore:
        expressions = ignore.split(',')
        for ex in expressions:
            files = [file for file in files if not ex in file]
    print('These files were found: ', files)
    table = []
    for i, file in enumerate(files):
        try:
            with open(file,'rb') as f:
                data = pickle.load(f)
            data['samples'] = samples_utils.unique_samples(data['samples'][-60:])
            parts = file.split('/')[-1].split('_samples.')[0].split('_')
            assert(commit==parts[0])
            ini = parts[1]
            i_blocks = parts.index('blocks')
            participant = '_'.join(parts[2:i_blocks])
            epoch = '_'.join([parts[i_blocks+1],parts[i_blocks+2]])
            if 'markov' in parts:
                model = 'markov'
            elif 'hmm8' in parts:
                model = 'hmm8'
            else:
                model = 'iHMM'

            for parameter in ['tau0','mu','sigma','jll','k']:
                if parameter == 'k':
                  p = [len(sample['phi']) for sample in data['samples']]
                else:
                  p = [sample[parameter] for sample in data['samples']]
                row = {'participant': participant,
                       'epoch': epoch,
                       'ini': ini,
                       'parameter': parameter,
                       'model': model,
                       'mean': np.mean(p),
                       'sd': np.std(p),
                       'q05': np.quantile(p,0.05),
                       'q95': np.quantile(p,0.95),
                       'n_eff': len(np.unique(p)),
                       'n_samples': len(p)}
                table.append(row)
            print('Done ({0:2.0f}%): '.format(i/len(files)*100), file)
        except Exception as e:
            print('Failed: ', file)
            print(e)
            traceback.print_exc()
    df = pd.DataFrame(table)
    if len(df)>0:
        output = pj(OUTPUT_DIR,commit,'sample_stats.csv')
        df.to_csv(output,float_format="%g", index=False)
        print('Sample stats saved to:', output)

def main(args):
    if args.parameter_settings:
        parameter_settings('parameter_settings.csv')
    if args.sample_stats:
        sample_stats(args.sample_stats, args.ignore)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Compute statistics for diagnostics.')
    parser.add_argument('--parameter_settings', action='store_true')
    parser.add_argument('--sample_stats', type=str, default=None)
    parser.add_argument('--ignore', type=str, default=None)
    args = parser.parse_args()
    main(args)