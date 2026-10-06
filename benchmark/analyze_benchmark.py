"""Summaries and performance profiles of the extended benchmark.

The method is used, as in the paper, within a multistart: from each official
starting point of an instance, ten perturbed copies are solved and the best
value is kept. A test problem is therefore an (instance, start) pair. A
variant solves it when the best value f* over its ten runs that stopped at
Step 5 is within a relative tolerance of the best value reached by any
variant; its cost is the total number of function evaluations (or CPU time)
of the ten runs. Performance profiles of Dolan and More are built on these
costs.

Usage:
    python analyze_benchmark.py [--variants FO H GN QN H1] [--tag main]
"""

import argparse
import collections
import csv
import math
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
RES = os.path.join(HERE, "results")

FTOL = 1.0e-3     # relative tolerance on the best f* for a problem to count as solved
LABELS = {"FO": "First-order ($B=0$)",
          "H": r"Shifted Hessian, $M=10^4$",
          "GN": r"Gauss--Newton, $M=10^4$",
          "QN": r"Shared BFGS, $M=10^4$",
          "H1": r"Shifted Hessian, $M=1$"}
STYLES = {"FO": ("k", "-"), "H": ("C0", (0, (1, 1))), "GN": ("C1", "-."),
          "QN": ("C2", ":"), "H1": ("C3", "--")}


def load():
    rows = list(csv.DictReader(open(os.path.join(RES, "benchmark.csv"))))
    for r in rows:
        for k in ("o", "start", "trial", "status", "subfail", "iterations", "nevals", "ncorrect"):
            r[k] = int(r[k])
        for k in ("cpu", "fstar"):
            r[k] = float(r[k])
    return rows


def problems(rows, variants):
    """Map (instance, start) -> {variant: list of runs}, complete problems only."""
    P = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        if r["variant"] in variants:
            P[(r["instance"], r["start"])][r["variant"]].append(r)
    ntr = max(len(v) for d in P.values() for v in d.values())
    return {k: d for k, d in P.items() if all(len(d[s]) == ntr for s in variants)}


def best_f(runs):
    vals = [r["fstar"] for r in runs if r["status"] == 0]
    return min(vals) if vals else math.inf


def outcome(P, variants):
    """For each problem and variant: (solved, total #fcnt, total CPU, identified)."""
    out = {}
    for key, d in P.items():
        fb = min(best_f(d[s]) for s in variants)
        out[key] = {}
        for s in variants:
            bf = best_f(d[s])
            ok = bf <= fb + FTOL * max(abs(fb), 1e-12)
            # outliers identified by the best run of the multistart
            ident = False
            if math.isfinite(bf):
                rb = min((r for r in d[s] if r["status"] == 0), key=lambda r: r["fstar"])
                ident = rb["ncorrect"] == rb["o"]
            out[key][s] = (ok, sum(r["nevals"] for r in d[s]),
                           sum(r["cpu"] for r in d[s]), ident, bf)
    return out


def ratios(O, variants, k):
    R = []
    for res in O.values():
        cost = [res[s][k] if res[s][0] else math.inf for s in variants]
        best = min(cost)
        if math.isfinite(best):
            R.append([c / best for c in cost])
    return np.array(R)


def plot_profiles(O, variants, tag):
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.3))
    for ax, k, title in zip(axes, (1, 2), ("Function evaluations", "CPU time")):
        R = ratios(O, variants, k)
        finite = R[np.isfinite(R)]
        tmax = max(2.0, float(np.max(finite)) * 1.05)
        tau = np.logspace(0, math.log10(tmax), 600)
        for j, s in enumerate(variants):
            rho = [(R[:, j] <= t).mean() for t in tau]
            c, ls = STYLES[s]
            ax.plot(tau, rho, color=c, linestyle=ls, label=LABELS[s], lw=1.6)
        ax.set_xscale("log", base=2)
        ax.set_xlim(1, tmax)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel(r"$\tau$")
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"$\rho_s(\tau)$")
    axes[1].legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    out = os.path.join(RES, "profiles_%s.pdf" % tag)
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


def summary(rows, O, variants, tag):
    meta = {}
    for f in os.listdir(os.path.join(HERE, "instances")):
        with open(os.path.join(HERE, "instances", f)) as fh:
            fh.readline()
            mid, n, m, k = map(int, fh.readline().split())
            meta[f[:-4]] = (n, m, k)
    lines = []
    head = "%-16s %2s %3s %2s" % ("problem", "n", "m", "o")
    for s in variants:
        head += " | %-3s %9s %8s" % (s, "f*", "#fcnt")
    lines.append(head)
    for key in sorted(O, key=lambda k: (meta[k[0]][1], k)):
        n, m, k = meta[key[0]]
        line = "%-16s %2d %3d %2d" % ("%s/s%d" % key, n, m, k)
        for s in variants:
            ok, fc, cpu, ident, bf = O[key][s]
            mark = ("*" if ok else " ") + ("i" if ident else " ")
            line += " | %s %9.3e %8d" % (mark, bf, fc)
        lines.append(line)
    lines.append("(* solved: best f* within %.0e of the best variant; i: all outliers identified)" % FTOL)
    stats = collections.defaultdict(lambda: [0, 0, 0])
    for res in O.values():
        for s in variants:
            stats[s][0] += res[s][0]
            stats[s][1] += res[s][3]
            stats[s][2] += res[s][1]
    for s in variants:
        st = [r for r in rows if r["variant"] == s and (r["instance"], r["start"]) in O]
        maxit = sum(r["status"] != 0 for r in st)
        sub = sum(r["subfail"] > 0 for r in st)
        lines.append("%-3s solved %2d/%d  identified %2d/%d  total #fcnt %9d  runs at maxit %3d  runs with inner-loop cap %3d"
                     % (s, stats[s][0], len(O), stats[s][1], len(O), stats[s][2], maxit, sub))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(RES, "summary_%s.txt" % tag), "w").write(txt + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="*", default=["FO", "H", "GN", "QN", "H1"])
    ap.add_argument("--tag", default="main")
    a = ap.parse_args()
    rows = load()
    P = problems(rows, a.variants)
    O = outcome(P, a.variants)
    print("%d test problems (instance, start)" % len(O))
    summary(rows, O, a.variants, a.tag)
    plot_profiles(O, a.variants, a.tag)


if __name__ == "__main__":
    main()
