import matlab.engine
import matlab
import numpy as np
print('Starting Matlab engine.')
eng = matlab.engine.start_matlab()
print('Finished starting Matlab engine.')

def rng(seed):
    eng.rng(seed)

def gamrnd(shape, scale, size=None):
    if np.size(shape) > 1:
        shape = shape.tolist()
        shape = matlab.double(shape)
    if np.size(scale) > 1:
        scale = scale.tolist()
        scale = matlab.double(scale)

    if type(shape) == np.float64:
        shape = np.float(shape)
    if type(scale) == np.float64:
        scale = np.float(scale) 

    if (size == None):
        return(np.array(eng.gamrnd(shape, scale)))
    else:
        return(np.array(eng.gamrnd(shape, scale, 1, size))[0])

def betarnd(a, b, size = None):
    if (np.size(a) > 1):
        a = a.tolist()
        a = matlab.double(a)
    if (np.size(b) > 1):
        b = b.tolist()
        b = matlab.double(b)

    if (type(a) == np.float64) or (type(a) == np.int64):
        a = np.float(a)
    if (type(b) == np.float64) or (type(b) == np.int64):
        b = np.float(b) 

    if (size == None) :
        return(np.array(eng.betarnd(a, b)))
    else:
        return(np.array(eng.betarnd(a, b, 1, matlab.double([size])))[0])

def binornd(n, p):
    if (np.size(n) > 1):
        n = n.tolist()
        n = matlab.double(n)
    if (np.size(p) > 1):
        p = p.tolist()
        p = matlab.double(p)

    return(np.array(eng.binornd(n, p)))

def rand():
    return(eng.rand())

def uniform(a, b):
    return(eng.rand()*(b-a)+a)