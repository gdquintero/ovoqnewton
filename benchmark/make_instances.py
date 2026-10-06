"""Build the OVO test instances used in the extended numerical study.

For each NIST StRD nonlinear regression problem, the observed data are
contaminated with a prescribed fraction of outliers, placed at random
observations and displaced by a fixed fraction of the range of the fitted
curve (alternating above and below it), so that they are clearly separated
from the noise of the clean data. The two NIST starting points and the
certified parameters are kept. The polynomial and Bard problems of the paper
are written in the same format with their original data.

Instance file format (whitespace separated):
    name
    model_id  n  m  nout
    lower bounds (n values)
    upper bounds (n values)
    start 1 (n values)
    start 2 (n values)
    reference parameters (n values)
    m lines: t  y
    nout lines: index (1-based) of each outlier
"""

import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
NIST_DIR = os.path.join(HERE, "nist")
OUT_DIR = os.path.join(HERE, "instances")
DATA_DIR = os.path.join(HERE, "..", "data")

LEVELS = [0.05, 0.10]
SEED = 20261005

# name: (model id, parameters whose sign must be kept in the box)
NIST = {
    "Misra1a":  (3, "all"),
    "Chwirut2": (4, "all"),
    "Rat43":    (5, "all"),
    "MGH17":    (6, [4, 5]),
    "Lanczos3": (7, "all"),
    "Thurber":  (8, [5, 6, 7]),
    "Kirby2":   (9, [4, 5]),
    "Bennett5": (10, [2, 3]),
    "ENSO":     (11, [4, 7]),
    "Gauss1":   (12, "all"),
}


def model(mid, t, b):
    """Plain-float version of the models of models.f90 (used to place outliers)."""
    if mid == 1:
        return b[0] + b[1] * t + b[2] * t ** 2 + b[3] * t ** 3
    if mid == 2:
        u, v = t, 16.0 - t
        return b[0] + u / (v * b[1] + min(u, v) * b[2])
    if mid == 3:
        return b[0] * (1.0 - math.exp(-b[1] * t))
    if mid == 4:
        return math.exp(-b[0] * t) / (b[1] + b[2] * t)
    if mid == 5:
        return b[0] / (1.0 + math.exp(b[1] - b[2] * t)) ** (1.0 / b[3])
    if mid == 6:
        return b[0] + b[1] * math.exp(-t * b[3]) + b[2] * math.exp(-t * b[4])
    if mid == 7:
        return b[0] * math.exp(-b[1] * t) + b[2] * math.exp(-b[3] * t) + b[4] * math.exp(-b[5] * t)
    if mid == 8:
        return (b[0] + b[1] * t + b[2] * t ** 2 + b[3] * t ** 3) / (
            1.0 + b[4] * t + b[5] * t ** 2 + b[6] * t ** 3)
    if mid == 9:
        return (b[0] + b[1] * t + b[2] * t ** 2) / (1.0 + b[3] * t + b[4] * t ** 2)
    if mid == 10:
        return b[0] * (b[1] + t) ** (-1.0 / b[2])
    if mid == 11:
        w = 2.0 * math.pi * t
        return (b[0] + b[1] * math.cos(w / 12) + b[2] * math.sin(w / 12)
                + b[4] * math.cos(w / b[3]) + b[5] * math.sin(w / b[3])
                + b[7] * math.cos(w / b[6]) + b[8] * math.sin(w / b[6]))
    if mid == 12:
        return (b[0] * math.exp(-b[1] * t) + b[2] * math.exp(-(t - b[3]) ** 2 / b[4] ** 2)
                + b[5] * math.exp(-(t - b[6]) ** 2 / b[7] ** 2))
    raise ValueError(mid)


def read_nist(name):
    """Return (start1, start2, certified, certified RSS, t, y) from a NIST .dat file."""
    lines = open(os.path.join(NIST_DIR, name + ".dat")).read().splitlines()
    s1, s2, cert = [], [], []
    rss = None
    data_start = None
    for k, line in enumerate(lines):
        tok = line.split()
        if len(tok) >= 5 and tok[0].startswith("b") and tok[1] == "=":
            s1.append(float(tok[2]))
            s2.append(float(tok[3]))
            cert.append(float(tok[4]))
        if line.strip().startswith("Residual Sum of Squares"):
            rss = float(tok[-1])
        if line.strip().startswith("Data:") and " y" in line:
            data_start = k + 1
    t, y = [], []
    for line in lines[data_start:]:
        tok = line.split()
        if len(tok) == 2:
            y.append(float(tok[0]))
            t.append(float(tok[1]))
    return s1, s2, cert, rss, t, y


def box(s1, s2, cert, keep_sign):
    """Box containing the starting points and the certified solution.

    Each interval is the hull of the three values, enlarged by half the
    largest magnitude. For parameters that must keep their sign (rates,
    widths, periods and denominators), the bound on the side of zero is
    replaced by one tenth of the value closest to zero.
    """
    lo, up = [], []
    for j, vals in enumerate(zip(s1, s2, cert)):
        a, b = min(vals), max(vals)
        pad = 0.5 * max(abs(v) for v in vals)
        lj, uj = a - pad, b + pad
        if keep_sign == "all" or (j + 1) in keep_sign:
            if a > 0:
                lj = 0.1 * a
            elif b < 0:
                uj = 0.1 * b
        lo.append(lj)
        up.append(uj)
    return lo, up


def contaminate(mid, cert, rss, t, y, level, rng):
    """Displace a fraction `level` of the observations (alternating signs)."""
    m = len(t)
    n = len(cert)
    k = max(2, int(round(level * m)))
    fit = [model(mid, ti, cert) for ti in t]
    span = max(fit) - min(fit)
    sres = math.sqrt(rss / max(1, m - n))
    idx = sorted(rng.sample(range(m), k))
    yc = list(y)
    for j, i in enumerate(idx):
        sign = 1.0 if j % 2 == 0 else -1.0
        yc[i] = y[i] + sign * (rng.uniform(0.2, 0.4) * span + 10.0 * sres)
    return yc, [i + 1 for i in idx]


def write_instance(fname, name, mid, lo, up, s1, s2, ref, t, y, outl):
    fmt = lambda v: " ".join("%.16e" % x for x in v)
    with open(os.path.join(OUT_DIR, fname), "w") as f:
        f.write(name + "\n")
        f.write("%d %d %d %d\n" % (mid, len(ref), len(t), len(outl)))
        for v in (lo, up, s1, s2, ref):
            f.write(fmt(v) + "\n")
        for ti, yi in zip(t, y):
            f.write("%.16e %.16e\n" % (ti, yi))
        for i in outl:
            f.write("%d\n" % i)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = random.Random(SEED)
    summary = []

    for name, (mid, keep) in NIST.items():
        s1, s2, cert, rss, t, y = read_nist(name)
        # Sanity check: the certified parameters must reproduce the certified RSS
        rss_chk = sum((model(mid, ti, cert) - yi) ** 2 for ti, yi in zip(t, y))
        assert abs(rss_chk - rss) <= 1e-6 * max(1.0, rss), (name, rss_chk, rss)
        lo, up = box(s1, s2, cert, keep)
        for level in LEVELS:
            yc, outl = contaminate(mid, cert, rss, t, y, level, rng)
            tag = "%s_%02d" % (name, int(round(100 * level)))
            write_instance(tag + ".txt", tag, mid, lo, up, s1, s2, cert, t, yc, outl)
            summary.append((tag, len(cert), len(t), len(outl)))

    # Polynomial of Section 4.1: 46 observations, the 10 with y = 10 are outliers
    rows = [tuple(map(float, l.split())) for l in open(os.path.join(DATA_DIR, "andreani.txt")).read().split("\n")[1:] if l.strip()]
    t = [r[0] for r in rows]
    y = [r[1] for r in rows]
    outl = [i + 1 for i, yi in enumerate(y) if yi == 10.0]
    x0 = [6.4602, 2.7072, -7.5418, 2.1604]
    write_instance("Poly_ref.txt", "Poly_ref", 1, [-10.0] * 4, [10.0] * 4, x0, x0, x0, t, y, outl)
    summary.append(("Poly_ref", 4, len(t), len(outl)))

    # Bard function of Section 4.2: 15 MGH observations plus 4 outliers
    rows = [tuple(map(float, l.split())) for l in open(os.path.join(DATA_DIR, "bard.txt")).read().split("\n")[1:] if l.strip()]
    t = [r[0] for r in rows]
    y = [r[1] for r in rows]
    outl = [16, 17, 18, 19]
    x0 = [1.0, 1.0, 1.0]
    write_instance("Bard_ref.txt", "Bard_ref", 2, [-10.0, 0.1, 0.1], [10.0] * 3, x0, x0, x0, t, y, outl)
    summary.append(("Bard_ref", 3, len(t), len(outl)))

    for tag, n, m, k in summary:
        print("%-12s n=%2d  m=%3d  outliers=%2d" % (tag, n, m, k))


if __name__ == "__main__":
    main()
