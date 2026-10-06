"""Run the extended benchmark of Algorithm 2.1 in parallel.

Each job runs one curvature variant on one instance from one NIST starting
point, with a multistart of perturbed copies of that point. All variants of
the same (instance, start) use the same random seed, hence the same starting
points. Jobs are pinned to distinct performance cores so that CPU times are
comparable across runs.

Usage:
    python3 run_benchmark.py [--trials 10] [--workers 8] [--only NAME ...]
                             [--variants FO H ...] [--cores 0 2 ...]
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "src", "ovo_bench")
INST = os.path.join(HERE, "instances")
RAW = os.path.join(HERE, "results", "raw")

# One logical CPU per performance core of the i7-13700KF (cores 0-7)
PCORES = [0, 2, 4, 6, 8, 10, 12, 14]

# Curvature variants: label -> driver arguments
VARIANTS = {
    "FO": ["B=0"],
    "H":  ["B=H", "M=1e4"],
    "GN": ["B=GN", "M=1e4"],
    "QN": ["B=QN", "M=1e4"],
    "H1": ["B=H", "M=1"],
    # Sensitivity to M (Comment 7): rescaling and eigenvalue clipping
    "H0.1": ["B=H", "M=0.1"],
    "H10": ["B=H", "M=10"],
    "H100": ["B=H", "M=100"],
    "Hc1": ["B=H", "M=1", "mmode=clip"],
    "Hc10": ["B=H", "M=10", "mmode=clip"],
    "GN1": ["B=GN", "M=1"],
}
MAIN = ["FO", "H", "GN", "QN", "H1"]


def jobs(trials, only, variants):
    for fname in sorted(os.listdir(INST)):
        tag = fname[:-4]
        if only and not any(tag.startswith(o) for o in only):
            continue
        ref = tag.endswith("_ref")
        starts = [1] if ref else [1, 2]
        for start in starts:
            seed = 123456 + 1000 * start
            common = ["inst=" + os.path.join(INST, fname), "start=%d" % start,
                      "ntrials=%d" % trials, "seed=%d" % seed, "maxit=10000", "project=1"]
            # The polynomial and Bard keep the perturbation of the paper
            common += ["pscale=abs", "pert=0.5"] if ref else ["pscale=rel", "pert=0.25"]
            for label in variants:
                args = VARIANTS[label]
                out = os.path.join(RAW, "%s_s%d_%s.csv" % (tag, start, label))
                yield label, out, common + args + ["out=" + out]


def run(job, cpu):
    label, out, args = job
    if os.path.exists(out):
        os.remove(out)
    cmd = ["taskset", "-c", str(cpu), EXE] + args
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    return out, res.returncode, res.stderr[-300:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=10)
    ap.add_argument("--workers", type=int, default=len(PCORES))
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--variants", nargs="*", default=MAIN)
    ap.add_argument("--cores", nargs="*", type=int, default=PCORES)
    a = ap.parse_args()
    cores = a.cores[: a.workers]

    os.makedirs(RAW, exist_ok=True)
    todo = list(jobs(a.trials, a.only, a.variants))
    # Longest jobs first improves load balance
    heavy = ("Gauss1", "Rat43", "ENSO", "Kirby2", "Bennett5")
    todo.sort(key=lambda j: not any(h in j[1] for h in heavy))

    free = list(cores)
    done = 0
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        pending = {}
        it = iter(todo)
        for job in it:
            cpu = free.pop()
            pending[pool.submit(run, job, cpu)] = cpu
            if not free:
                break
        while pending:
            fut = next(as_completed(pending))
            cpu = pending.pop(fut)
            out, rc, err = fut.result()
            done += 1
            status = "ok" if rc == 0 else "FAILED (%d) %s" % (rc, err)
            print("[%3d/%3d] %s %s" % (done, len(todo), os.path.basename(out), status), flush=True)
            nxt = next(it, None)
            if nxt is not None:
                pending[pool.submit(run, nxt, cpu)] = cpu

    # Merge all raw files into a single table
    merged = os.path.join(HERE, "results", "benchmark.csv")
    header = ("instance,o,B,M,mmode,tau,delta,sigmin,start,trial,status,subfail,"
              "iterations,cpu,fstar,nevals,ncorrect,variant\n")
    with open(merged, "w") as f:
        f.write(header)
        for name in sorted(os.listdir(RAW)):
            label = name[:-4].split("_")[-1]
            for line in open(os.path.join(RAW, name)):
                f.write(line.rstrip("\n") + "," + label + "\n")
    print("merged into", merged)


if __name__ == "__main__":
    sys.exit(main())
