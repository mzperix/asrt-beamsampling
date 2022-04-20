import pandas as pd
import glob
from os.path import join as pj
import numpy as np
import pickle
from definitions import ROOT_DIR, DATA_DIR, CONFIG_DIR, TEMPLATE_DIR
import configparser
import argparse

def get_experiment_data(experiment_name):
    if experiment_name.lower() == 'twm':
        data = pd.read_csv(pj(DATA_DIR, 'ASRT_TWM_merged.txt'),
                           delimiter ='\t',
                           encoding = 'iso8859_2')

    if experiment_name.lower() == 'asrt_180':
        data = pd.read_csv(pj(DATA_DIR, 'ASRT_180.csv'))
        data = data.rename(columns={'participant': 'Subject',
                                    'session': 'Session',
                                    'block': 'Block',
                                    'stimulus': 'event',
                                    'rt': 'firstRT',
                                    'correct_response': 'firstACC',
                                    'key': 'firstRESP', # CHECK IF THIS IS CORRECT
                                    'trial_type': 'TrialType',
                                    })
        #print(data.head())
        data['Trial'] = np.tile(range(1,86), len(data)//85)
        data['Block'] = data['Block']+(data['Session']-1)*15

    if experiment_name.lower() == 'elarasztas':
        data = pd.read_csv(pj(DATA_DIR,'elarasztas_dataset.csv'))
        # # Session 1-8
        # data = pd.read_csv(pj(DATA_DIR,'random_learn_iszt.txt'), encoding = 'ISO-8859-2', delimiter='\t')
        # # Session 9
        # df = pd.read_csv(pj(DATA_DIR,'isztivel.txt'), encoding = 'ISO-8859-2', delimiter='\t')
        # # Change IDs in session 9
        # ids = pd.read_csv(pj(DATA_DIR,'elarasztas_ID_azonositok.csv'), delimiter='\t')
        # df['Subject'] = df[['Subject']].replace(to_replace=ids.set_index('ID_session9')[['ID_session10']].to_dict()['ID_session10'])
        # data = data.append(df)
        # # Session 10
        # data = data.append(pd.read_csv(pj(DATA_DIR,'ti8zeshez02.txt'), encoding = 'ISO-8859-2', delimiter='\t'))
        # data['Trial'] = np.tile(range(1,86), len(data)//85)
        # data['Block'] = data['Block']+(data['Session']-1)*25

    #if experiment_name.lower() == 'napexp_teszt':
    #    data = pd.read_csv(pj(data_dir,'asrt_teszt.txt'),
    #                                delimiter ='\t',
    #                                encoding = 'iso8859_2')
    #if experiment_name.lower() == 'napexp_reteszt':
    #    data = pd.read_csv(pj(data_dir,'asrt_reteszt.txt'),
    #                       delimiter ='\t',
    #                       encoding = 'iso8859_2')

    if experiment_name.lower() == 'implicit_erp':
        data = pd.read_csv(pj(DATA_DIR, 'alldata_behav_implicit_ASRT_n32.csv'), delimiter=';')
        data = data.rename(columns={
            'ID': 'Subject',
            'scenario': 'Session',
            'block': 'Block',
            'proba': 'Trial',
            'position': 'event',
            'RT': 'firstRT',
            'talalat': 'firstACC',
            'resp': 'firstRESP', # CHECK IF THIS IS CORRECT
            'triplet': 'TT',
            'trialtype': 'TrialType',
        })
        data['TrialType'] = data.TrialType.replace({1:'P',2:'R'})
        data['Subject'] = data.Subject.str.slice(stop=2)

    data = filter_columns(data)
    return(data)

def filter_columns(data):
    return(data[['Subject','Session','Block', 'firstRT','event','firstACC',
                                                  'TrialType', 'firstRESP','Trial']])

def get_unique_values(experiment_name):
    data = get_experiment_data(experiment_name)
    print('--- Participants ---')
    print(','.join(map(str,data['Subject'].unique())))
    print('\n---   Sessions   ---')
    print(','.join(map(str,data['Session'].unique())))
    print('\n---    Blocks    ---')
    print(','.join(map(str,data['Block'].unique())))
    print('\n---    Trials    ---')
    print(','.join(map(str,data['Trial'].unique())))
    print('\n---Session&Block ---')
    print(data.groupby(['Session','Block']).size())
    print('\n--- Trial Types  ---')
    print(','.join(map(str,data['TrialType'].unique())))
    print('\n---   Stimuli    ---')
    print(','.join(map(str,data['event'].unique())))

def get_participant_data(data, participant):
    return(filter_columns(data.loc[data.Subject == participant]))

def get_ini_list():
    return([filename.split('/')[-1][:-4] for filename in glob.glob(pj(CONFIG_DIR,'*.ini'))])

def get_ini(filename):
    if '_markov' in filename:
        filename = filename.split('_')[0]
    elif '_' in filename:
        filename = filename.split('_')[-1]
    if '.ini' not in filename:
        filename = filename+'.ini'
    config = configparser.ConfigParser()
    config.read(pj(CONFIG_DIR,filename))
    ini = dict(**config)
    for key,value in ini.items():
        ini[key] = dict(**value)
    return(ini)

def get_hypers(filename):
    ini = get_ini(filename)
    hypers = ini['HYPERS']
    return(hypers)

def data_filename(participant, block):
    return(participant+'_blocks_'+block+'.pkl')

def import_data(**kwargs):
    """
    Imports datafiles from DATA_DIR (given in definitions).
    USAGE:
        - filename: name of datafile
        OR
        - participant: participant (experiment ID and participant ID)
        - block: blocks (e.g. 16_20)
    """
    if 'filename' in kwargs:
        with open(pj(DATA_DIR, kwargs['filename']),'rb') as file:
            return(pickle.load(file))
    if ('participant' in kwargs) and ('block' in kwargs):
        with open(pj(DATA_DIR, data_filename(kwargs['participant'], kwargs['block'])), 'rb') as file:
            return(pickle.load(file))
    raise NotImplemented('Use filename or participant and block.')

def export_blocks(experiment_name, dataset, participant, blocks_from, blocks_to, **kwargs):
    data = get_participant_data(dataset, participant)
    data = data.loc[data['Block']<=blocks_to]
    data = data.loc[data['Block']>=blocks_from]

    output = dict(Y = np.array(data.event.tolist())-1,
                  rt = np.array(data.firstRT.tolist()),
                  correct_response = np.array(data.firstACC.tolist()),
                  trial_type = np.array(data.TrialType.tolist()),
                  block = np.array(data.Block.tolist()),
                  first_response = np.array(data.firstRESP.tolist()),
                  trial = np.array(data.Trial.tolist()),
                  )

    # Extra options
    #output = filter_by_rt(output, **kwargs)

    if len(output['Y'])>0:
        filename = experiment_name+'_'+str(participant)+'_blocks_'+str(blocks_from)+'_'+str(blocks_to)+'.pkl'
        with open(pj(DATA_DIR,filename), 'wb') as file:
            pickle.dump(output, file)
        print('Done for participant', participant, 'blocks: ', str(blocks_from),'-',str(blocks_to))
    else:
        print('Participant',participant,'not found.')

    return(filename)

def export_blocks_with_input():
    experiment = input('Experiment name:')
    participants = input('Participant IDs (with comma separation):')
    blocks_from = input('Blocks from (with comma separation):')
    blocks_to = input('Blocks to (with comma separation, same number as blocks_from):')
    try:
        # GENERATE CORRESPONDING INPUT DATA
        print('Loading experiment data.')
        data = get_experiment_data(experiment)
        if participants == 'all':
            participants = data.Subject.unique()
        else:
            participants = [int(p) for p in participants.split(',')]
        if ',' in blocks_from:
            blocks_from = [int(b) for b in blocks_from.split(',')]
            blocks_to = [int(b) for b in blocks_to.split(',')]

        print('Generating input files.')
        input_files = []
        for p in participants:
            for block_from, block_to in zip(blocks_from, blocks_to):
                input_files.append(export_blocks(experiment, data, p, block_from, block_to))
        print(input_files)

    except Exception as e:
        print('An unexpected error occurred.')
        print(e)

def filter_by_rt(data, rt_min=0, rt_max=0):
    if 'filter' not in data:
        data['filter'] = np.ones(np.size(data['Y']))
    data['filter'] = data['filter'] * (data['rt']>=rt_min)
    if rt_max > 0:
        data['filter'] = data['filter'] * (data['rt']<=rt_max)
    return(data)

def print_stats(data):
    print('---- ALL ----')
    print('RT mean, sd: ', np.mean(data['rt']), np.std(data['rt']))
    print('RT min, max: ', np.min(data['rt']), np.max(data['rt']))
    print('---- CORRECT ----')
    print('RT mean, sd: ', np.mean(data['rt'][data['correct_response']==1]), np.std(data['rt'][data['correct_response']==1]))
    print('RT min, max: ', np.min(data['rt'][data['correct_response']==1]), np.max(data['rt'][data['correct_response']==1]))
    print('---- INCORRECT ----')
    print('RT mean, sd: ', np.mean(data['rt'][data['correct_response']==0]), np.std(data['rt'][data['correct_response']==0]))
    print('RT min, max: ', np.min(data['rt'][data['correct_response']==0]), np.max(data['rt'][data['correct_response']==0]))

def print_data(data):
    print(pd.DataFrame(data))

def get_template(template_name):
    config = configparser.ConfigParser()
    config.read(pj(TEMPLATE_DIR,template_name+'.template'))
    template = dict(**config)
    for key,value in template.items():
        template[key] = dict(**value)
    return(template)

def main():
    parser = argparse.ArgumentParser(description='Data handling.')
    parser.add_argument('-f', metavar = 'file', type=str,
                        help='Filename')
    parser.add_argument('--stats', action='store_true',
                        help='Print out stats')
    parser.add_argument('--data', action='store_true',
                        help='Print out data')
    parser.add_argument('-u', type=str, help='Print list of unique column values of a given experiment.')
    parser.add_argument('--convert', action='store_true',
                        help='Convert raw experiment data into input for cogtom.py')
    args = parser.parse_args()

    if args.f is not None:
        with open(args.f, 'rb') as file:
            data = pickle.load(file)
        if args.stats:
            print_stats(data)
        if args.data:
            print_data(data)
    if args.u:
        get_unique_values(args.u)

    if args.convert:
        export_blocks_with_input()


if __name__=='__main__':
    main()
