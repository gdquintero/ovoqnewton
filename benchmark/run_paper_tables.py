"""Rerun the experiments of Sections 4.1-4.4 of the paper with ovo_bench.

Same settings as the submitted version (100 perturbed starting points, seed
123456, M = 1, delta = sigma_min = 0.1, epsilon = 1e-4), with the stopping
test of Step 5 and the starting points projected onto the box, as stated in
the paper. Results are written to results/paper/*.csv.

Usage:
    python3 run_paper_tables.py [--workers 8] [--skip-sens]
"""

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "src", "ovo_bench")
INST = os.path.join(HERE, "instances")
OUT = os.path.join(HERE, "results", "paper")
PCORES = [0, 2, 4, 6, 8, 10, 14, 12]


def job_list(skip_sens):
    base = ["ntrials=100", "seed=123456", "pscale=abs", "pert=0.5", "project=1",
            "M=1", "maxit=10000"]
    jobs = []
    # Tables 2 and 4 (polynomial, o = 0..12) and Tables 3 and 4 (Bard, o = 0..6)
    for inst, omax in (("Poly_ref", 12), ("Bard_ref", 6)):
        for B in ("0", "H"):
            for o in range(omax + 1):
                out = os.path.join(OUT, "%s_B%s_o%02d.csv" % (inst, B, o))
                jobs.append((out, base + ["inst=%s/%s.txt" % (INST, inst), "B=" + B,
                                          "o=%d" % o, "out=" + out]))
    # Table 5: sensitivity of the Bard experiment to delta and sigma_min
    if not skip_sens:
        for delta in ("1e-1", "1e-2", "1e-3", "1e-4"):
            for sigmin in ("1e-1", "1e-2", "1e-3"):
                for B in ("0", "H"):
                    for o in range(1, 7):
                        out = os.path.join(OUT, "sens_d%s_s%s_B%s_o%d.csv" % (delta, sigmin, B, o))
                        jobs.append((out, base + ["inst=%s/Bard_ref.txt" % INST, "B=" + B,
                                                  "o=%d" % o, "delta=" + delta,
                                                  "sigmin=" + sigmin, "out=" + out]))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=len(PCORES))
    ap.add_argument("--skip-sens", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    jobs = job_list(a.skip_sens)

    def run(k):
        out, args = jobs[k]
        if os.path.exists(out):
            os.remove(out)
        cpu = PCORES[k % a.workers]
        subprocess.run(["taskset", "-c", str(cpu), EXE] + args, check=True,
                       capture_output=True)
        return out

    # Jobs are distributed round-robin, so each worker keeps its own core
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures = []
        for w in range(a.workers):
            idx = list(range(w, len(jobs), a.workers))
            futures.append(pool.submit(lambda ids: [run(k) for k in ids], idx))
        for f in futures:
            f.result()
    print("done:", len(jobs), "jobs")


if __name__ == "__main__":
    main()
