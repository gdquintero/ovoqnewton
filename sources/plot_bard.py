import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import models
import os

cwd = os.getcwd()
parent = os.path.abspath(os.path.join(cwd, os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4, size_img * 4.8]
plt.rcParams['axes.unicode_minus'] = False
plt.rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['cmr10', 'STIXGeneral', 'DejaVu Serif']

df_data = pd.read_table(parent+"/data/bard.txt", delimiter=" ", header=None, skiprows=1, skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_bard.txt", delimiter=" ", header=None, skiprows=0, skipinitialspace=True)
df_outliers = pd.read_table(parent+"/output/outliers_bard.txt", delimiter=" ", header=None, skiprows=0, skipinitialspace=True)

# The model is evaluated over the design range [1,15]; the data are not sorted
# by abscissa (the outliers were appended), so the curve range is set explicitly.
t = np.linspace(1, 15, 1000)
noutliers = df_outliers[0].values[0]
outliers = np.empty((2, noutliers))

for i in range(noutliers):
    outliers[0, i] = df_data[0].values[df_outliers[0].values[i+1]-1]
    outliers[1, i] = df_data[1].values[df_outliers[0].values[i+1]-1]

plt.plot(df_data[0].values, df_data[1].values, "ko", ms=2)
plt.plot(t, models.bard(t, *df_sol.values[0]), lw=1)
plt.plot(outliers[0], outliers[1], 'ro', mfc='none', ms=6, mew=0.6)
plt.tick_params(axis='both', direction='in')
plt.xlabel(r'$t$')
plt.ylabel(r'$y$')
plt.xticks(np.arange(1, 16, 2))
plt.yticks(np.arange(-2, 5, 1))
plt.xlim(0.5, 15.5)
plt.ylim(-2.5, 5)
plt.savefig(parent+"/images/bard_fitting.pdf", bbox_inches="tight")
plt.show()
