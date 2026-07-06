import numpy as np

def andreani(t,x1,x2,x3,x4):
    return x1 + (x2 * t) + (x3 * (t**2)) + (x4 * (t**3))

def meyer(t,x1,x2,x3):
    return x1 * np.exp(x2 / (t + x3))

def bard(t,x1,x2,x3):
    # t is the observation index i; u=i, v=16-i, w=min(u,v)
    u = t
    v = 16.0 - t
    w = np.minimum(u,v)
    return x1 + u / (v*x2 + w*x3)

def osborne(t,x1,x2,x3,x4,x5):
    # Osborne 1 (MGH #17): sum of two decaying exponentials, n=5, m=33.
    return x1 + x2*np.exp(-t*x4) + x3*np.exp(-t*x5)