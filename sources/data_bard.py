import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import models
import os

cwd = os.getcwd()
parent = os.path.abspath(os.path.join(cwd, os.pardir))

seed = 42
rng = np.random.default_rng(seed)

# 15 clean Bard observations; the abscissa is the observation index i (u=i,
# v=16-i, w=min). New outlier points reuse the same design via their abscissa.
t_clean = np.arange(1, 16, dtype=float)
y_clean = np.array([
    0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39,
    0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39
])

# Least-squares fit to the clean data, used only to place the outliers relative
# to the model curve.
popt, _ = curve_fit(models.bard, t_clean, y_clean, p0=[1.0, 1.0, 1.0], maxfev=20000)

# Append 4 outliers at random abscissas, two above and two below the curve, with
# random magnitudes clearly separated from the clean-data noise (reproducible via
# the fixed seed).
n_out = 4
t_out = rng.uniform(2.0, 14.0, n_out)
signs = np.array([1.0, 1.0, -1.0, -1.0])          # two above, two below
mag = rng.uniform(1.5, 2.5, n_out)
y_out = models.bard(t_out, *popt) + signs * mag

t = np.concatenate([t_clean, t_out])
y = np.concatenate([y_clean, y_out])
m = t.size                                          # 19

with open(parent+"/data/bard.txt", "w") as f:
    f.write("%i\n" % m)
    for i in range(m):
        f.write("%f %f\n" % (t[i], y[i]))

tt = np.linspace(1, 15, 1000)
plt.plot(t_clean, y_clean, "ko", ms=3)
plt.plot(t_out, y_out, "rs", ms=5)
plt.plot(tt, models.bard(tt, *popt), lw=1)
plt.show()
