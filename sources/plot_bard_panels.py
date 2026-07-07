import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import models
import os

cwd = os.getcwd()
parent = os.path.abspath(os.path.join(cwd, os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['axes.unicode_minus'] = False
# Use mathtext with Computer Modern instead of usetex: it matches the paper's
# font and, unlike the usetex PDF backend, embeds the minus sign correctly.
plt.rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['cmr10', 'STIXGeneral', 'DejaVu Serif']

fig, (ax1, ax2) = plt.subplots(1, 2,
    figsize=[size_img * 6.4 * 2.0, size_img * 4.8])

# --- Panel (a): fitted model at o = 4 ---
df_data = pd.read_table(parent+"/data/bard.txt", delimiter=" ", header=None,
                        skiprows=1, skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_bard.txt", delimiter=" ",
                       header=None, skiprows=0, skipinitialspace=True)
df_outliers = pd.read_table(parent+"/output/outliers_bard.txt", delimiter=" ",
                            header=None, skiprows=0, skipinitialspace=True)

t = np.linspace(1, 15, 1000)
noutliers = df_outliers[0].values[0]
outliers = np.empty((2, noutliers))
for i in range(noutliers):
    outliers[0, i] = df_data[0].values[df_outliers[0].values[i+1]-1]
    outliers[1, i] = df_data[1].values[df_outliers[0].values[i+1]-1]

ax1.plot(df_data[0].values, df_data[1].values, "ko", ms=2)
ax1.plot(t, models.bard(t, *df_sol.values[0]), lw=1)
ax1.plot(outliers[0], outliers[1], 'ro', mfc='none', ms=6, mew=0.6)
ax1.tick_params(axis='both', direction='in')
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$y$')
ax1.set_xticks(np.arange(1, 16, 2))
ax1.set_yticks(np.arange(-2, 5, 1))
ax1.set_xlim(0.5, 15.5)
ax1.set_ylim(-2.5, 5)
ax1.set_title(r'(a)', y=-0.32)

# --- Panel (b): best f(x*) as a function of o (semilog) ---
df_drop = pd.read_table(parent+"/output/drop_bard.txt", delimiter=" ",
                        header=None, skipinitialspace=True)
o = df_drop[0].values
fx = df_drop[1].values

ax2.semilogy(o, fx, 'k.-', lw=1, ms=6)
ax2.tick_params(axis='both', which='both', direction='in')
ax2.set_xlabel(r'Number of outliers $o$')
ax2.set_ylabel(r'$f(x^*)$')
ax2.set_xticks(np.arange(0, 7, 1))
ax2.set_ylim(3e-5, 5.0)
ax2.set_yticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
ax2.set_yticklabels([r'$10^{0}$', r'$10^{-1}$', r'$10^{-2}$',
                     r'$10^{-3}$', r'$10^{-4}$'])
ax2.set_title(r'(b)', y=-0.32)

fig.tight_layout()
fig.savefig(parent+"/images/bard_panels.pdf", bbox_inches="tight")
plt.show()
