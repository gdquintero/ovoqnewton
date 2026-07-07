import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

# o and best f(x*) over the 100 starting points (the two variants coincide).
df = pd.read_table(parent+"/output/drop_bard.txt", delimiter=" ",
                   header=None, skipinitialspace=True)
o  = df[0].values
fx = df[1].values

plt.semilogy(o, fx, 'k.-', lw=1, ms=6)
plt.tick_params(axis='both', which='both', direction='in')
plt.xlabel(r'Number of outliers $o$')
plt.ylabel(r'$f(x^*)$')
plt.xticks(np.arange(0, 7, 1))
plt.ylim(3e-5, 5.0)
plt.yticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4],
           [r'$10^{0}$', r'$10^{-1}$', r'$10^{-2}$', r'$10^{-3}$', r'$10^{-4}$'])
plt.savefig(parent+"/images/bard_drop.pdf", bbox_inches="tight")
plt.show()
