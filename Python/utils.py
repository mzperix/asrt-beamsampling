import os
import sys
import subprocess
import numpy as np
import itertools
import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool
import time
#from multiprocessing import Pool
from definitions import ROOT_DIR
#import git

def find_file(path):
#    for _dir in os.environ['PATH'].split(':'):
    for _dir in sys.path:
        if _dir == '':
            if os.path.exists(path):
                return(os.path.abspath(path))
        else:
            if os.path.exists(_dir+'/'+path):
                return(os.path.abspath(_dir+'/'+path))
        print(_dir+'/'+path)
    return('File not found.')

def get_git_hash(path=ROOT_DIR):
    try:
        output = subprocess.check_output(['cat',os.path.join(path,'.git/refs/heads/recovery')])
        return(output.decode().split('\n')[0][:7])
    except:
        return('NOGTHSH')


def print_log_matrix(matrix):
    symbols = [' ', '.', 'o', 'O']
    l_matrix = np.floor(np.log10(matrix))
    for row in l_matrix:
        print('[', end='')
        for e in row:
            if e<=-4:
                print(symbols[0], end='')
            else:
                print(symbols[int(e)], end='')
        print(']')

def arg_cube(*args):
    new_args = []
    N = 1
    for arg in args:
        if type(arg) == list:
            new_args.append(arg)
            N *= len(arg)
        else:
            new_args.append([arg])
    arglist = itertools.product(*new_args)
    return(arglist, N)

def kwarg_cube(**kwargs):
    args = []
    keys = []
    new_kwargs = []
    for k,v in kwargs.items():
        keys.append(k)
        args.append(v)
    cube_args, N = arg_cube(*args)
    for arg in cube_args:
        kwarg = dict()
        for k,v in zip(keys,arg):
            kwarg[k] = v
        new_kwargs.append(kwarg)
    return(new_kwargs, N)

def pd_multiplicator(func):
    def wrapper(*args, **kwargs):
        i = 0
        progress = 0
        def progress_dot():
            while progress < (i*10)//(N*M):
                print('.',end='')
                progress += 1

        def arg_extend(a):
            #print('.',end='')
            return(func(*a[0],**a[1]))

        if 'parallel' in kwargs:
            parallel = kwargs.pop('parallel')
        else:
            parallel = 1
        df = pd.DataFrame()
        arglist, N = arg_cube(*args)
        kwarglist, M = kwarg_cube(**kwargs)
        arguments = []
        for arg in arglist:
            for kwarg in kwarglist:
                arguments.append((arg,kwarg))
        if parallel > 1:
            pool = Pool(nodes=parallel)
            results = []
            start_time = time.time()
            for i, res in enumerate(pool.imap(arg_extend, arguments)):
                elapsed_time = time.time() - start_time
                remaining_time = elapsed_time*(len(arguments)-i-1)/(i+1)
                sys.stderr.write('\r{0:2.4%}  {1} elapsed. Remaining: {2}'\
                    .format(i/len(arguments),
                        time.strftime("%H:%M:%S", time.gmtime(elapsed_time)),
                        time.strftime("%H:%M:%S", time.gmtime(remaining_time))))
                results.append(res)
            print('')
            #results = pool.imap_unordered(arg_extend, arguments)
            df = pd.concat(results)
            #pool.terminate()
        else:
            print('[',end='')
            for arg, kwarg in arguments:
                i += 1
                df = df.append(func(*arg,**kwarg))
                while progress < (i*10)//(N*M):
                    print('.',end='')
                    progress += 1
            print(']')
        return(df)
    return wrapper

def ton(d, key):
    if d and (key in d):
        return d[key]
    else:
        return None
