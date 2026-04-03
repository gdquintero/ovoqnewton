import numpy as np

l = lambda t, a, b, c: (a * t - c) * np.exp(-b * t) + c

F = lambda t, a, b, c: 1.0 - np.exp((a/b) * t * np.exp(-b * t) + (1.0/b) * ((a/b) - c) * \
    (np.exp(-b * t) - 1.0) - c * t)

def osborne2(t,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11):
    y = x1 * np.exp(-t * x5) + x2 * np.exp(-x6 * (t - x9)**2) + \
        x3 * np.exp(-x7 * (t - x10)**2) + x4 * np.exp(-x8 * (t - x11)**2)
    
    return y

def andreani(t,x1,x2,x3,x4):
    return x1 + (x2 * t) + (x3 * (t**2)) + (x4 * (t**3))