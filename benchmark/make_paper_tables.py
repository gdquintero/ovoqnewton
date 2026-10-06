"""LaTeX rows of the tables of Sections 4.1-4.4, from results/paper/*.csv.

Prints the body rows of the polynomial and Bard tables, the CPU-time table and
the sensitivity table, in the format used in the manuscript.
"""

import csv
import glob
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, "results", "paper")


def load(name):
    return list(csv.reader(open(os.path.join(P, name))))


def block(rows):
    """best f, #it and #fcnt of the best run, mean/median #it and #fcnt, mean/median CPU."""
    best = min(rows, key=lambda r: float(r[14]))
    it = [int(r[12]) for r in rows]
    fc = [int(r[15]) for r in rows]
    cpu = [1e3 * float(r[13]) for r in rows]
    return (float(best[14]), int(best[12]), int(best[15]),
            st.mean(it), st.median(it), st.mean(fc), st.median(fc),
            st.mean(cpu), st.median(cpu))


def sweep_table(inst, omax, label):
    print("%% %s" % label)
    tot = {B: [0] * 9 for B in "0H"}
    for o in range(omax + 1):
        cells = []
        for B in "0H":
            b = block(load("%s_B%s_o%02d.csv" % (inst, B, o)))
            for k in range(1, 9):
                tot[B][k] += b[k]
            cells.append("%.3e & %d & %d & %.0f & %.1f & %.0f & %.1f" % b[:7])
        print("%-2d & %s & %s \\\\" % (o, cells[0], cells[1]))
    print("\\hline")
    print("$\\Sigma$ & & %d & %d & %.0f & %.1f & %.0f & %.1f & & %d & %d & %.0f & %.1f & %.0f & %.1f \\\\"
          % (tuple(tot["0"][1:7]) + tuple(tot["H"][1:7])))
    return tot


def cpu_sum(inst, os_):
    out = {}
    for B in "0H":
        mean = sum(block(load("%s_B%s_o%02d.csv" % (inst, B, o)))[7] for o in os_)
        med = sum(block(load("%s_B%s_o%02d.csv" % (inst, B, o)))[8] for o in os_)
        out[B] = (mean, med)
    return out


def main():
    sweep_table("Poly_ref", 12, "Table: polynomial")
    sweep_table("Bard_ref", 6, "Table: Bard")

    print("% Table: CPU time (ms)")
    for lab, inst, os_ in (("Polynomial, $o=0,\\dots,12$", "Poly_ref", range(13)),
                           ("Bard, $o=0,\\dots,6$", "Bard_ref", range(7)),
                           ("Bard, $o=1,\\dots,6$", "Bard_ref", range(1, 7))):
        c = cpu_sum(inst, os_)
        print("%s\n  & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$ & $%.2f$ & $%.2f$ \\\\"
              % (lab, c["0"][0], c["0"][1], c["H"][0], c["H"][1],
                 c["H"][0] / c["0"][0], c["H"][1] / c["0"][1]))

    print("% Table: sensitivity (sums over o = 1..6 of the mean counts)")
    wins_fc = wins_it = 0
    for d in ("1e-1", "1e-2", "1e-3", "1e-4"):
        for s in ("1e-1", "1e-2", "1e-3"):
            cnt = {}
            for B in "0H":
                it = fc = 0.0
                for o in range(1, 7):
                    b = block(load("sens_d%s_s%s_B%s_o%d.csv" % (d, s, B, o)))
                    it += b[3]
                    fc += b[5]
                cnt[B] = (it, fc)
            f3 = block(load("sens_d%s_s%s_B0_o3.csv" % (d, s)))[0]
            f4 = block(load("sens_d%s_s%s_B0_o4.csv" % (d, s)))[0]
            wins_fc += cnt["H"][1] < cnt["0"][1]
            wins_it += cnt["H"][0] < cnt["0"][0]
            tex = lambda v: "$10^{%d}$" % int(v.split("e")[1])
            print("%s & %s & %.3e & %.3e & $%.0f\\times$ & %.0f & %.0f & %.0f & %.0f \\\\"
                  % (tex(d), tex(s), f3, f4, f3 / f4, cnt["0"][0], cnt["0"][1], cnt["H"][0], cnt["H"][1]))
    print("%% second-order lower mean #fcnt in %d/12, lower mean #it in %d/12" % (wins_fc, wins_it))


if __name__ == "__main__":
    main()
