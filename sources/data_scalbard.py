import sys
import numpy as np
import models
import os

# Generate synthetic data for the Bard scalability experiment.
# Usage: python data_scalbard.py <m>
if len(sys.argv) < 2:
    sys.exit("Usage: python data_scalbard.py <m>   (e.g. 100, 1000, 10000)")

m = int(sys.argv[1])

cwd = os.getcwd()
parent = os.path.abspath(os.path.join(cwd, os.pardir))

# Fixed seed so the dataset is reproducible for each m and identical across the
# first-order and quasi-Newton runs.
rng = np.random.default_rng(42)

# Ground-truth parameters: the classical Bard minimizer. With it the clean curve
# reproduces the m=15 Bard data of Section 4.2 (y in ~[0.13, 4.40]).
xtrue = np.array([0.08241, 1.133, 2.344])

# Abscissas in [1,15], so that v = 16 - t and w = min(t, 16 - t) as in Section 4.2.
t = np.linspace(1.0, 15.0, m)
yclean = models.bard(t, *xtrue)

# Each observation is an outlier with probability 0.1; if so, it lies above the
# curve with probability 0.8 (up to 8) and below with probability 0.2 (down to -3).
is_out = rng.random(m) < 0.1
above = rng.random(m) < 0.8

y = yclean + rng.uniform(-0.1, 0.1, m)                       # clean: y + U[-0.1,0.1]
y_above = rng.uniform(yclean, 8.0)                           # outlier above: U[y,8]
y_below = rng.uniform(-3.0, yclean)                          # outlier below: U[-3,y]
y = np.where(is_out & above, y_above, y)
y = np.where(is_out & ~above, y_below, y)

obar = int(is_out.sum())

with open(parent+"/data/scalbard.txt", "w") as f:
    f.write("%i %i\n" % (m, obar))
    for i in range(m):
        f.write("%f %f\n" % (t[i], y[i]))

print("data_scalbard: m=%i, generated outliers obar=%i" % (m, obar))
