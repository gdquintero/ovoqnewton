"""Analysis of the active-set study (run_active.py).

For each instance, delta and variant: mean |I_delta(x^k)| over the runs, CPU
time per subproblem (total CPU over the number of trial points, each of which
requires one subproblem and one evaluation of f), and the ratio of the
second-order variants to the first-order method. Prints LaTeX rows and writes
the figure results/active_cost.pdf (CPU time per subproblem versus the mean
size of the active set).
"""

import csv
import glob
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11})
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["cmr10", "STIXGeneral", "DejaVu Serif"]
plt.rcParams["axes.formatter.use_mathtext"] = True

HERE = os.path.dirname(os.path.abspath(__file__))
ACT = os.path.join(HERE, "results", "active")
INST = [("Bennett5_10", "Bennett5", 3), ("Kirby2_10", "Kirby2", 5),
        ("ENSO_10", "ENSO", 9), ("Gauss1_10", "Gauss1", 8)]
DELTAS = ["1e-1", "1", "1e1", "1e2", "1e3"]
VARS = ["FO", "H1", "GN1"]


def stats(inst, d, v):
    rows = list(csv.reader(open(os.path.join(ACT, "%s_d%s_%s.csv" % (inst, d, v)))))
    it = np.array([int(r[12]) for r in rows])
    cpu = np.array([float(r[13]) for r in rows])
    fc = np.array([int(r[15]) for r in rows])
    msz = np.array([float(r[17]) for r in rows])
    # Each evaluation after the first corresponds to one subproblem solve
    return {"m": float(np.mean(msz)), "mspi": 1e3 * cpu.sum() / max(1, it.sum()),
            "msps": 1e3 * cpu.sum() / max(1, (fc - 1).sum()),
            "it": float(np.mean(it)), "cpu": float(cpu.sum())}


def main():
    S = {(i, d, v): stats(i, d, v) for i, _, _ in INST for d in DELTAS for v in VARS}
    print("% instance & delta & |I| & ms/subproblem FO & H & GN & H/FO & GN/FO")
    for inst, lab, n in INST:
        for d in DELTAS:
            fo, h, gn = (S[(inst, d, v)] for v in VARS)
            tex = "$10^{%d}$" % int(float(d) and round(np.log10(float(d))))
            print("%s & %s & %.1f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\"
                  % (lab if d == DELTAS[0] else "", tex, fo["m"], fo["msps"], h["msps"],
                     gn["msps"], h["msps"] / fo["msps"], gn["msps"] / fo["msps"]))

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    marks = {"Bennett5_10": "o", "Kirby2_10": "s", "ENSO_10": "^", "Gauss1_10": "D"}
    cols = {"FO": "k", "H1": "C3", "GN1": "C1"}
    names = {"FO": "First-order", "H1": r"Shifted Hessian, $M=1$", "GN1": r"Gauss--Newton, $M=1$"}
    for v in VARS:
        for inst, lab, n in INST:
            xs = [S[(inst, d, "FO")]["m"] for d in DELTAS]
            ys = [S[(inst, d, v)]["msps"] for d in DELTAS]
            ax.loglog(xs, ys, color=cols[v], marker=marks[inst], ms=4, lw=0.8,
                      mfc="none" if v != "FO" else cols[v])
    for v in VARS:
        ax.plot([], [], color=cols[v], lw=1, label=names[v])
    for inst, lab, n in INST:
        ax.plot([], [], "k", marker=marks[inst], ls="none", ms=4, mfc="none",
                label="%s ($n=%d$)" % (lab, n))
    ax.set_xlabel(r"mean size of $I_\delta(x^k)$")
    ax.set_ylabel("CPU time per subproblem (ms)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.5, ncol=2, loc="upper left", frameon=False)
    fig.tight_layout()
    out = os.path.join(HERE, "results", "active_cost.pdf")
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
