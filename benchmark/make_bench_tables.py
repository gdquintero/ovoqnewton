"""LaTeX rows for the tables of Sections 4.5 and 4.6 (extended benchmark).

For each set of variants, a test problem (instance, start) counts as solved
by a variant when its best value is within 1e-3 (relative) of the best value
reached by any variant of the same set.
"""

import analyze_benchmark as A

MAIN = ["FO", "H1", "H", "GN", "QN"]
MSET = ["FO", "H0.1", "H1", "H10", "H100", "H", "Hc1", "Hc10", "HA1", "HA10", "HA100",
        "GN1", "GN", "QN"]
DESC = {
    "FO":   ("First-order", "--", "--"),
    "H0.1": ("Shifted Hessian", "$0.1$", "rescale"),
    "H1":   ("Shifted Hessian", "$1$", "rescale"),
    "H10":  ("Shifted Hessian", "$10$", "rescale"),
    "H100": ("Shifted Hessian", "$10^2$", "rescale"),
    "H":    ("Shifted Hessian", "$10^4$", "rescale"),
    "Hc1":  ("Shifted Hessian", "$1$", "clip"),
    "Hc10": ("Shifted Hessian", "$10$", "clip"),
    "HA1":  ("Shifted Hessian", "$\\sigma_{k,j}$", "adaptive"),
    "HA10": ("Shifted Hessian", "$10\\,\\sigma_{k,j}$", "adaptive"),
    "HA100": ("Shifted Hessian", "$10^2\\sigma_{k,j}$", "adaptive"),
    "GN1":  ("Gauss--Newton", "$1$", "rescale"),
    "GN":   ("Gauss--Newton", "$10^4$", "rescale"),
    "QN":   ("Shared BFGS", "$10^4$", "rescale"),
}


def table(rows, variants, title):
    P = A.problems(rows, variants)
    O = A.outcome(P, variants)
    print("%% %s: %d test problems" % (title, len(O)))
    for s in variants:
        runs = [r for r in rows if r["variant"] == s and (r["instance"], r["start"]) in O]
        solved = sum(O[k][s][0] for k in O)
        ident = sum(O[k][s][3] for k in O)
        fc = sum(O[k][s][1] for k in O)
        cpu = sum(O[k][s][2] for k in O)
        maxit = sum(r["status"] == 1 for r in runs)
        tlim = sum(r["status"] >= 2 for r in runs)
        cap = sum(r["subfail"] > 0 for r in runs)
        name, M, mode = DESC[s]
        print("%s & %s & %s & %d & %d & %.2f & %.1f & %d & %d \\\\"
              % (name, M, mode, solved, ident, fc / 1e6, cpu / 60.0, maxit + tlim, cap))


def main():
    rows = A.load()
    table(rows, MAIN, "main benchmark")
    table(rows, MSET, "sensitivity to M")


if __name__ == "__main__":
    main()
