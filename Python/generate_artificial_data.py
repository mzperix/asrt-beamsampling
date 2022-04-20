import sys
sys.path.append('Python')
import pickle
import numpy as np
import carpenter_williams_model as cwm
import beam_sampling as bs
from ideal_observer import IdealObserver
import configparser
import argparse
import pandas
import itertools

ASRT = 'ASRT'
RANDOM = 'RANDOM'

def asrt_emission(state):
    if state % 2 == 0:
        return(int(state/2))
    else:
        return(np.random.choice(range(4)))

def generate_asrt_input(N, *args):
    T = N*8
    S0 = np.tile(np.array([0,1,2,3,4,5,6,7]),N)
    Y = np.array([asrt_emission(s) for s in S0])
    trial_type = np.tile(['P','R'],T//2)
    return(Y, trial_type)

def generate_random_sequence(N, D):
    T = N*8
    Y = np.random.randint(D, size=T)
    np.repeat('R', T)
    return(Y, trial_type)

def generate_sequence(sequence_type, *args):
    if sequence_type==ASRT:
        Y, trial_type = generate_asrt_input(*args)
    if sequence_type==RANDOM:
        Y, trial_type = generate_random_sequence(*args)
    return(Y, trial_type)

def sample_internal_model(N, steps, sequence_type=ASRT, *args):
    print('Sampling internal model:  N =',N, ' steps =',steps)
    Y, trial_type = generate_sequence(sequence_type, N, *args)
    print('Y:', Y)
    K0 = 40
    S00 = np.array([np.random.choice(range(K0)) for n in range(N*8)])
    hypers = dict(alpha0 = 3.84, gamma = 1.3,
                  alpha0_a = 4.0, alpha0_b = 1.0,
                  gamma_a = 4.0, gamma_b = 1.0,
                  H = np.array([[1.0,1.0,1.0,1.0]]))
    numb, nums, numi = steps, 1, 1
    samples, stats = bs.sample_beam_ihmm(Y.copy(), hypers.copy(), numb, nums, numi, S00.copy())
    return(samples[-1]['Pi'], samples[-1]['Phi'])

def generate_config():
    cfgfile = open("params.ini",'w')

    Config = configparser.ConfigParser()
    # add the settings to the structure of the file, and lets write it out...
    Config.add_section('INTERNAL_MODEL')
    Config.set('INTERNAL_MODEL','seed', '1')
    Config.set('INTERNAL_MODEL','N', '20')
    Config.set('INTERNAL_MODEL','steps', '700')

    Config.add_section('RT_SAMPLING')
    Config.set('RT_SAMPLING','seed', '1')
    Config.set('RT_SAMPLING','N', '30')
    Config.set('RT_SAMPLING','tau0', '4.0')
    Config.set('RT_SAMPLING','mu', '0.012')
    Config.set('RT_SAMPLING','sigma', '0.0005')
    Config.write(cfgfile)
    cfgfile.close()

def generate_csv(output, cube=False):
    param_names = ['tau0', 'mu', 'sigma', 'internal_model_n', 'internal_model_steps', 'internal_model_seed', 'rt_n', 'rt_n_blocks', 'rt_seed',]
    df = pandas.DataFrame(columns=param_names)
    if cube:
        print('Create an n-dimensional grid for parameter sweep.')
        print('CAUTION: the value of internal_model_n and rt_n is multiplied by 8 (full sequence) to get the number of trials.')
        print('         rt_n means number of cycles within one block')
        params = dict()
        for name in param_names:
            lower = np.float32(input('Lower bound of '+name+':'))
            upper = np.float32(input('Upper bound of '+name+':'))
            if name[-4:] == 'seed':
                n = int(upper-lower+1)
            else:
                n = int(input('Resolution of '+name+':'))
            params[name]=np.linspace(lower, upper, n)
        df = pandas.DataFrame(list(itertools.product(*params.values())), columns=[*params.keys()])
    df.to_csv(output, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate_config", action="store_true",
                        help="Generate config file for parameters of artificial data")
    parser.add_argument("--generate_csv", action="store_true",
                        help="Generate csv file for parameters of artificial data")
    parser.add_argument("--cube", action="store_true",
                        help="Ask parameter ranges for generate_csv and create row for all possible parameter settings")
    parser.add_argument("--rt_csv", action="store_true", help="Generate rt")
    parser.add_argument("--rt_ini", action="store_true", help="Generate rt")
    parser.add_argument("--input", type = str, help = "Input config file")
    parser.add_argument("--output", type = str, help = "Output filename (pickle)")
    parser.add_argument("--seq", type=str, default='ASRT',
        help = "Type of sequence for training and sequence generation for RT production.")
    parser.add_argument("-d", type=int, default=4,
        help = "Number of alternative observations.")
    args = parser.parse_args()
    config = configparser.ConfigParser()

    if args.generate_config:
        generate_config()

    if args.generate_csv:
        generate_csv(args.output, args.cube)

    elif args.rt_ini:
        print('Currently not working.')
        pass
        input_file = args.input
        output_file = args.output

        config = configparser.ConfigParser()
        config.read(input_file)

        rt_params = config['RT_SAMPLING']
        np.random.seed(int(config['INTERNAL_MODEL']['seed']))
        Pi, Phi = sample_internal_model(int(config['INTERNAL_MODEL']['N']), int(config['INTERNAL_MODEL']['steps']))
        observer = IdealObserver(tau0 = np.float32(rt_params['tau0']),
                                 mu = np.float32(rt_params['mu']),
                                 sigma = np.float32(rt_params['sigma']),
                                 Pi = Pi,
                                 Phi = Phi)

        np.random.seed(int(rt_params['seed']))
        Y = generate_asrt_input(int(rt_params['N']))
        s_hat, rt_log_pred_prob, rt = observer.generate_rt(Y)

        with open(output_file, 'wb') as file:
            pickle.dump(output_data, file)

    elif args.rt_csv:
        sequence_type = args.seq
        D = args.d
        model_cache = dict()
        output_file = args.output
        table = pandas.read_csv(args.input, header = 0,
                               dtype = {'internal_model_seed': np.int32,
                                        'internal_model_n': np.int32,
                                        'rt_seed': np.int32})
        ID = 0
        for r in table.iterrows():
            row = r[1]
            print('Generating rt with following parameters: ')
            print(row)
            model_id = (sequence_type+'_'+
                        str(row['internal_model_n'])+'_'+
                        str(row['internal_model_seed'])+'_'+
                        str(row['internal_model_steps']))
            if model_id in model_cache:
                Pi = model_cache[model_id]['Pi']
                Phi = model_cache[model_id]['Phi']
            else:
                np.random.seed(int(row['internal_model_seed']))
                Pi, Phi = sample_internal_model(
                    int(row['internal_model_n']),
                    int(row['internal_model_steps']),
                    sequence_type,
                    D)
                model_cache[model_id] = dict(internal_model_n = row['internal_model_n'],
                                             internal_model_seed = row['internal_model_seed'],
                                             internal_model_steps = row['internal_model_steps'],
                                             Pi = Pi,
                                             Phi = Phi)
                print(model_cache)
            observer = IdealObserver(tau0 = row['tau0'],
                                     mu = row['mu'],
                                     sigma = row['sigma'],
                                     Pi = Pi,
                                     Phi = Phi)

            np.random.seed(int(row['rt_seed']))
            n_blocks = int(row['rt_n_blocks'])
            rt_n = int(row['rt_n'])
            block_length = rt_n * 8
            Y, trial_type = generate_sequence(sequence_type, rt_n*n_blocks, D)
            s_hat, log_pred_prob, rt = observer.generate_rt(Y)
            block = np.repeat(list(range(1,n_blocks+1)),block_length)
            trial = np.tile(list(range(1,block_length+1)),n_blocks)
            while len(block)<len(Y):
                block.append((len(block)//block_length)+1)
                trial.append((len(trial)%block_length)+1)
            block = np.int32(block)
            trial = np.int32(trial)

            output_data = dict(rt = rt,
                               Y = Y,
                               block = block,
                               trial = trial,
                               correct_response = np.repeat(1, len(trial)),
                               filters = np.repeat(True,len(trial)),
                               trial_type = trial_type
                               )
            model_data = dict(log_pred_prob = log_pred_prob,
                              Pi = Pi,
                              Phi = Phi,
                              tau0 = observer.tau0,
                              mu = observer.mu,
                              sigma = observer.sigma,
                              K = np.size(Phi,0))
            output_rt = (output_file+'_'+sequence_type+'_'+str(ID)+'_blocks_'+
                str(np.min(block))+'_'+str(np.max(block))+'.pkl')
            output_model = (output_file+'_'+sequence_type+'_'+str(ID)+'.model')
            ID += 1
            with open(output_rt, 'wb') as file:
                pickle.dump(output_data, file)
            with open(output_model, 'wb') as file:
                pickle.dump(model_data, file)
            print(f'RT stats: min-max {min(rt), max(rt)} mean, sd {np.mean(rt), np.std(rt)}')
            print('Saved RT results into ',output_rt)
            print('Saved model into      ',output_model)

if __name__ == "__main__":
    main()
