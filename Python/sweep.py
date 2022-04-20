# Creating ini files for parameter sweeping
import configparser
from os.path import join as pj

data_dir = '../Data'

base_config = configparser.ConfigParser()
base_config.read(pj(data_dir,'curling.ini'))

H = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.3, 1.8, 2.5]
alpha0 = [1.3]
gamma = [3.8]
epsilon_min = [0.02, 0.005]
epsilon_max = [0.4, 0.3, 0.2, 0.1]

filename = input('Ini filenames: ')

counter = 0
for h in H:
    for a in alpha0:
        for g in gamma:
            for emin in epsilon_min:
                for emax in epsilon_max:
                    base_config.set('HYPERS', 'H', str(h))
                    base_config.set('HYPERS', 'alpha0', str(a))
                    base_config.set('HYPERS', 'gamma', str(g))
                    base_config.set('HYPERS', 'epsilon_min', str(emin))
                    base_config.set('HYPERS', 'epsilon_max', str(emax))
                    with open(pj(data_dir, 'curling_'+filename+str(counter)+'.ini'),'w') as fileobject:
                        base_config.write(fileobject, space_around_delimiters=True)
                    counter += 1