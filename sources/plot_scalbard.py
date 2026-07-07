import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
parent = os.path.abspath(os.path.join(cwd, os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['axes.unicode_minus'] = False
plt.rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['cmr10', 'STIXGeneral', 'DejaVu Serif']

sizes = [100, 1000, 10000]
labels = ['(a)', '(b)', '(c)']

fig, axes = plt.subplots(1, 3,
    figsize=[size_img * 6.4 * 3.0, size_img * 4.8])

for ax, m, lab in zip(axes, sizes, labels):
    df_fo = pd.read_table(parent+f"/output/scalbard_fo_curve_{m}.txt",
                          delim_whitespace=True, header=None)
    df_qn = pd.read_table(parent+f"/output/scalbard_qn_curve_{m}.txt",
                          delim_whitespace=True, header=None)

    ax.semilogy(df_fo[0].values, df_fo[1].values, 'k.-', lw=1, ms=4,
                label=r'$B=0$')
    ax.semilogy(df_qn[0].values, df_qn[1].values, 'r.--', lw=1, ms=3,
                label=r'$M=1$')
    ax.tick_params(axis='both', which='both', direction='in')
    ax.set_xlabel(r'Number of discarded observations $o$')
    if m == sizes[0]:
        ax.set_ylabel(r'$f(x^*)$')
    ax.set_title(lab + r'  $m=%s$' % ('{:,}'.format(m)), y=-0.36)
    ax.legend(frameon=False, fontsize=9)

fig.tight_layout()
fig.savefig(parent+"/images/scalbard.pdf", bbox_inches="tight")
print("saved images/scalbard.pdf")
