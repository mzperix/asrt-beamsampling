import pickle
import glob
import sys
import numpy as np
import os
import time
import utils

def main():
    if len(sys.argv)<2:
        print('No expression given')
        return()
    else:
        if sys.argv[1] == 'latest':
            expr = utils.get_git_hash('..')[:7]
        else:
            expr = sys.argv[1]
        files = np.sort(glob.glob('Output/*'+expr+'*'))
        L = np.max([len(f) for f in files])
        # print(files)
        for f in files:
            try:
                with open(f,'rb') as ff:
                    data = pickle.load(ff)
                if not 'markov' in f:
                    K = [sample['K'] for sample in data['samples']]
                print(f, end='')
                for k in range(L-len(f)):
                    print(' ', end='')
                print('   Number of samples: {:5d}'.format(len(data['samples'])), end = '')
                if not 'markov' in f:
                    print('   K min, max, mean: {:3d} {:3d} {:5.2f}'.format(np.min(K), np.max(K), np.mean(K)), 
                          end = '')
                if np.floor((time.time()-os.path.getmtime(f))/60) < 5:
                    print('    R')
                else:
                    print('')
            except:
                print('Error with file: ',f)
if __name__=='__main__':
    main()
