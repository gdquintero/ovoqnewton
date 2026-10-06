"""Adaptive bound ||B_{k,j,i}|| <= c * sigma_{k,j} on the extended benchmark.

Same test problems, seeds and settings as run_benchmark.py, for the shifted
Hessian with c in {1, 10, 100}. Each trial runs as a separate process
(first=k reproduces the k-th starting point), under a wall-clock limit, so
that a subproblem on which Algencan does not return cannot block the study;
such a trial is recorded with status 3. Results are appended to results/raw/
with the labels HA1, HA10 and HA100.

Usage:
    python3 run_adaptive.py [--workers 8]
"""

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "src", "ovo_bench")
INST = os.path.join(HERE, "instances")
RAW = os.path.join(HERE, "results", "raw")
TMP = os.path.join(HERE, "results", "adaptive_trials")
PCORES = [0, 2, 4, 6, 8, 10, 12, 14]
VARIANTS = {"HA1": ["B=H", "cad=1"], "HA10": ["B=H", "cad=10"], "HA100": ["B=H", "cad=100"]}
NTRIALS = 10
WALL = 2400   # seconds of wall-clock time per trial (the CPU limit is 1800 s)


def trial_jobs():
    for fname in sorted(os.listdir(INST)):
        tag = fname[:-4]
        ref = tag.endswith("_ref")
        for start in ([1] if ref else [1, 2]):
            seed = 123456 + 1000 * start
            common = ["inst=" + os.path.join(INST, fname), "start=%d" % start,
                      "seed=%d" % seed, "maxit=10000", "project=1", "maxtime=1800"]
            common += ["pscale=abs", "pert=0.5"] if ref else ["pscale=rel", "pert=0.25"]
            for label, args in VARIANTS.items():
                for k in range(1, NTRIALS + 1):
                    out = os.path.join(TMP, "%s_s%d_%s_t%02d.csv" % (tag, start, label, k))
                    yield (tag, start, label, k, out,
                           common + args + ["ntrials=%d" % k, "first=%d" % k, "out=" + out])


def run(job, cpu):
    tag, start, label, k, out, args = job
    if os.path.exists(out):
        os.remove(out)
    try:
        subprocess.run(["taskset", "-c", str(cpu), EXE] + args, capture_output=True,
                       timeout=WALL)
    except subprocess.TimeoutExpired:
        # The subproblem solver did not return: record an unsuccessful trial
        with open(out, "w") as f:
            f.write("%s,0,H,0,scale,0,0,0,%d,%d,3,0,0,%.5E,1.0E+99,0,0,0,0\n"
                    % (tag, start, k, float(WALL)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=len(PCORES))
    a = ap.parse_args()
    os.makedirs(TMP, exist_ok=True)
    jobs = list(trial_jobs())
    heavy = ("Gauss1", "Rat43", "ENSO", "Kirby2", "Bennett5")
    jobs.sort(key=lambda j: not any(h in j[0] for h in heavy))

    def worker(w):
        for idx in range(w, len(jobs), a.workers):
            run(jobs[idx], PCORES[w])

    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        list(pool.map(worker, range(a.workers)))

    # Assemble one raw file per (instance, start, variant), as run_benchmark.py
    groups = {}
    for tag, start, label, k, out, _ in jobs:
        groups.setdefault((tag, start, label), []).append((k, out))
    for (tag, start, label), items in groups.items():
        with open(os.path.join(RAW, "%s_s%d_%s.csv" % (tag, start, label)), "w") as f:
            for k, out in sorted(items):
                line = open(out).readline().rstrip("\n")
                # keep the 17 columns of the benchmark format
                f.write(",".join(line.split(",")[:17]) + "\n")
    print("done:", len(jobs), "trials")


if __name__ == "__main__":
    main()
