import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import models
import os

cwd = os.getcwd()
parent = os.path.abspath(os.path.join(cwd, os.pardir))

seed = 42
rng = np.random.default_rng(seed)

# 33 clean Osborne 1 observations; the abscissa is t_i = 10*(i-1).
t_clean = 10.0 * np.arange(0, 33, dtype=float)
y_clean = np.array([
    0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881, 0.850, 0.818,
    0.784, 0.751, 0.718, 0.685, 0.658, 0.628, 0.603, 0.580, 0.558,
    0.538, 0.522, 0.506, 0.490, 0.478, 0.467, 0.457, 0.448, 0.438,
    0.431, 0.424, 0.420, 0.414, 0.411, 0.406
])

# Least-squares fit to the clean data, used only to place the outliers relative
# to the model curve.
popt, _ = curve_fit(models.osborne, t_clean, y_clean,
                    p0=[0.5, 1.5, -1.0, 0.01, 0.02], maxfev=40000)

# Append 4 outliers at random abscissas, two above and two below the curve, with
# random magnitudes clearly separated from the clean-data noise (reproducible via
# the fixed seed). y-range of the clean data is ~0.4-0.94, so offsets ~0.3-0.5
# make unambiguous outliers.
n_out = 4
t_out = rng.uniform(30.0, 290.0, n_out)
signs = np.array([1.0, 1.0, -1.0, -1.0])          # two above, two below
mag = rng.uniform(0.3, 0.5, n_out)
y_out = models.osborne(t_out, *popt) + signs * mag

t = np.concatenate([t_clean, t_out])
y = np.concatenate([y_clean, y_out])
m = t.size                                          # 37

with open(parent+"/data/osborne.txt", "w") as f:
    f.write("%i\n" % m)
    for i in range(m):
        f.write("%f %f\n" % (t[i], y[i]))

print("LS fit (clean):", popt)
print("outlier abscissas:", t_out)

tt = np.linspace(0, 320, 1000)
plt.plot(t_clean, y_clean, "ko", ms=3)
plt.plot(t_out, y_out, "rs", ms=5)
plt.plot(tt, models.osborne(tt, *popt), lw=1)
plt.show()
