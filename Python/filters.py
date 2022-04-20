import numpy as np
import data_handling as dh

def correct(data):
    return(np.array(data['correct_response'])==1)

def rt_min(rt, rt_min):
    return(np.array(rt)>=rt_min)

def rt_max(rt, rt_max):
    return(np.array(rt)<=rt_max)

def skip_start_at_block(data,n):
    return(np.array(data['trial'])>n)
    
def new_stream(data):
    flags = [1]
    current_block = data['block'][0]
    for block in data['block'][1:]:
        if block == current_block:
            flags.append(0)
        else:
            flags.append(1)
            current_block = block
    flags = np.int8(flags)
    return(flags)
   
def get_filter(data,ini):
    if isinstance(ini,str):
        config = dh.get_ini(ini)
    elif isinstance(ini,dict):
        config = ini
    else:
        config = ini
    #else:
    #    raise NotImplemented('Please use string or dict ini.')
    filters = np.ones(len(data['rt']), dtype = bool)
    if 'FILTERS' in config:
        rt = data['rt']
        if 'rt_min' in config['FILTERS']:
            filters *= rt_min(rt,int(config['FILTERS']['rt_min']))
        if 'rt_max' in config['FILTERS']:
            filters *= rt_max(rt,int(config['FILTERS']['rt_max']))
        if 'skip_start_at_block' in config['FILTERS']:
            filters *= skip_start_at_block(data, int(config['FILTERS']['skip_start_at_block']))
        if 'correct_response' in config['FILTERS']:
            if config['FILTERS']['correct_response'] == '1':
                filters *= correct(data)
    return(filters)

def filter(data, ini):
    filters = get_filter(data, ini)
    return(np.array(data['rt'])[filters])