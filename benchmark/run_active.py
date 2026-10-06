"""Cost per iteration as a function of the size of the active set (Comment 5).

The band half-width delta controls the size of I_delta(x^k). On the four
largest instances, delta is increased from 0.1 to 1000 and the first-order
method, the shifted Hessian (M = 1) and the Gauss-Newton matrix (M = 1) are
run from the same ten perturbed starting points. The driver reports the mean
and maximum |I_delta(x^k)| of each run. Results go to results/active/.

Usage:
    python3 run_active.py [--workers 8]
"""

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "src", "ovo_bench")
OUT = os.path.join(HERE, "results", "active")
PCORES = [0, 2, 4, 6, 8, 10, 12, 14]

INSTANCES = ["Bennett5_10", "Kirby2_10", "ENSO_10", "Gauss1_10"]
DELTAS = ["1e-1", "1", "1e1", "1e2", "1e3"]
VARIANTS = {"FO": ["B=0"], "H1": ["B=H", "M=1"], "GN1": ["B=GN", "M=1"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=len(PCORES))
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    for inst in INSTANCES:
        for d in DELTAS:
            for lab, args in VARIANTS.items():
                out = os.path.join(OUT, "%s_d%s_%s.csv" % (inst, d, lab))
                jobs.append((out, ["inst=%s/instances/%s.txt" % (HERE, inst), "start=1",
                                   "ntrials=10", "seed=124456", "pscale=rel", "pert=0.25",
                                   "project=1", "maxit=10000", "maxtime=600",
                                   "delta=" + d, "out=" + out] + args))

    def worker(w):
        for k in range(w, len(jobs), a.workers):
            out, args = jobs[k]
            if os.path.exists(out):
                os.remove(out)
            subprocess.run(["taskset", "-c", str(PCORES[w]), EXE] + args,
                           capture_output=True)

    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        list(pool.map(worker, range(a.workers)))
    print("done:", len(jobs), "jobs")


if __name__ == "__main__":
    main()
