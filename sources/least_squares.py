import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import models
import os

cwd = os.getcwd()
parent = os.path.abspath(os.path.join(cwd, os.pardir))

# Problem name is passed as a command-line argument, e.g.:
#   python least_squares.py bard
if len(sys.argv) < 2:
    sys.exit("Usage: python least_squares.py <problem>   (e.g. bard, andreani)")

problem = sys.argv[1]

# Standard starting points per problem (default: all ones).
init_points = {
    "bard":     np.array([1.0, 1.0, 1.0]),
    "andreani": np.array([1.0, 1.0, 1.0, 1.0]),
}

# Box constraints per problem, matching the Fortran driver. For Bard-type models
# the box x2,x3 in [0.1,10] keeps the denominator positive and excludes the
# spurious minimum (x2,x3 -> +-inf) into which an unconstrained fit diverges.
bounds = {
    "bard":     ([-10.0, 0.1, 0.1], [10.0, 10.0, 10.0]),
}

model = getattr(models, problem)
init = init_points.get(problem, np.ones(model.__code__.co_argcount - 1))

df = pd.read_table(parent+"/data/"+problem+".txt", delimiter=" ", header=None,
                   skiprows=1, skipinitialspace=True)

if problem in bounds:
    popt, pcov = curve_fit(model, df[0].values, df[1].values, p0=init,
                           bounds=bounds[problem], maxfev=20000)
else:
    popt, pcov = curve_fit(model, df[0].values, df[1].values, p0=init, maxfev=20000)

with open(parent+"/output/sol_ls_"+problem+".txt", "w") as f:
    for v in popt:
        f.write("%f\n" % v)

print("LS fit ("+problem+"):", popt)

t = np.linspace(df[0].values[0], df[0].values[-1], 1000)
plt.plot(df[0].values, df[1].values, "ko", ms=2)
plt.plot(t, model(t, *popt), lw=1)
plt.show()
