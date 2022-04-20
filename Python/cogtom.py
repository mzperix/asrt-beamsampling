import sys
sys.path.append('Notebooks/Carpenter-Williams')
#import pickle
import os
import argparse

def main():
    if 'PYTHON_START' in os.environ:
        os.chdir(os.environ['PYTHON_START'])

    if sys.argv[1] == 'markov':
        import markov_model
        if len(sys.argv)>3:
            markov_model.sample_markov_from_file(sys.argv[2], sys.argv[3], model_type=markov_model.LATER_MODEL)
        else:
            markov_model.sample_markov_from_file(sys.argv[2], model_type=markov_model.LATER_MODEL)

    elif sys.argv[1] == 'linear_markov':
        import markov_model
        if len(sys.argv)>3:
            markov_model.sample_markov_from_file(sys.argv[2], sys.argv[3], model_type=markov_model.LINEAR_MODEL)
        else:
            markov_model.sample_markov_from_file(sys.argv[2], model_type=markov_model.LINEAR_MODEL)

    elif sys.argv[1] == 'hmm':
        import hmm
        if len(sys.argv)<4:
            return('Not enough arguments passed. Need number of states')
        else:
            hmm.sample_hmm_from_file(sys.argv[2], sys.argv[3], K=int(sys.argv[4]))

    else:
        import beam_sampling_rt_s_hat as bs_rt
        if len(sys.argv) > 2:
            bs_rt.sample_beam_ihmm_from_file(sys.argv[1], sys.argv[2])
        else:
            bs_rt.sample_beam_ihmm_from_file(sys.argv[1])

if __name__ == "__main__":
    main()
