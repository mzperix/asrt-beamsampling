import pickle
import glob
from os.path import join as pj
import pandas as pd

DIR = "../Data/artificial_asrt"

for filename in glob.glob(pj(DIR,'*.model')):
    with open(filename, 'rb') as f:
        d = pickle.load(f)
    print(filename)

    participant_test = int(filename[:-6].split('_')[-1])
    output = filename[:-6]+'_probs.csv'
    pd.DataFrame({'participant_test': participant_test, 'log_pred_prob': d['log_pred_prob']}).to_csv(output)
