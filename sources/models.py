import numpy as np

def andreani(t,x1,x2,x3,x4):
    return x1 + (x2 * t) + (x3 * (t**2)) + (x4 * (t**3))

def meyer(t,x1,x2,x3):
    return x1 * np.exp(x2 / (t + x3)) 