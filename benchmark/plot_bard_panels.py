"""Figure of Section 4.2: Bard fit at o = 4 and best f(x*) versus o.

Panel (a) uses the fit of the first-order method at o = 4 (output/), which is
unchanged by the revision. Panel (b) plots the best f(x*) over the 100
starting points for both variants, from results/paper/.
"""

import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
size_img = 0.6

plt.rcParams.update({"font.size": 11})
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["cmr10", "STIXGeneral", "DejaVu Serif"]
plt.rcParams["axes.formatter.use_mathtext"] = True


def bard(t, x1, x2, x3):
    u, v = t, 16.0 - t
    return x1 + u / (v * x2 + np.minimum(u, v) * x3)


data = np.loadtxt(os.path.join(ROOT, "data", "bard.txt"), skiprows=1)
sol = np.loadtxt(os.path.join(ROOT, "output", "solution_bard.txt"))
out = np.loadtxt(os.path.join(ROOT, "output", "outliers_bard.txt"), dtype=int)
idx = out[1:] - 1

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[size_img * 6.4 * 2.0, size_img * 4.8])

t = np.linspace(1, 15, 1000)
ax1.plot(data[:, 0], data[:, 1], "ko", ms=2)
ax1.plot(t, bard(t, *sol), lw=1)
ax1.plot(data[idx, 0], data[idx, 1], "ro", mfc="none", ms=6, mew=0.6)
ax1.tick_params(axis="both", direction="in")
ax1.set_xlabel(r"$t$")
ax1.set_ylabel(r"$y$")
ax1.set_xticks(np.arange(1, 16, 2))
ax1.set_yticks(np.arange(-2, 5, 1))
ax1.set_xlim(0.5, 15.5)
ax1.set_ylim(-2.5, 5)
ax1.set_title(r"(a)", y=-0.32)

os_ = np.arange(7)
for B, style, lab in (("0", "k.-", r"First-order ($B=0$)"),
                      ("H", "C0s--", r"Second-order ($M=1$)")):
    best = []
    for o in os_:
        rows = csv.reader(open(os.path.join(HERE, "results", "paper", "Bard_ref_B%s_o%02d.csv" % (B, o))))
        best.append(min(float(r[14]) for r in rows))
    ax2.semilogy(os_, best, style, lw=1, ms=4 if B == "H" else 6, mfc="none" if B == "H" else None,
                 label=lab)
ax2.tick_params(axis="both", which="both", direction="in")
ax2.set_xlabel(r"Number of outliers $o$")
ax2.set_ylabel(r"$f(x^*)$")
ax2.set_xticks(os_)
ax2.set_ylim(3e-5, 5.0)
ax2.set_yticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
ax2.legend(fontsize=8, loc="lower left", frameon=False)
ax2.set_title(r"(b)", y=-0.32)

fig.tight_layout()
target = os.path.join(ROOT, "latex", "bard_panels.pdf")
fig.savefig(target, bbox_inches="tight")
print("wrote", target)
