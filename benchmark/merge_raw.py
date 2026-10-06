"""Merge results/raw/*.csv into results/benchmark.csv (as run_benchmark.py does)."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "results", "raw")
header = ("instance,o,B,M,mmode,tau,delta,sigmin,start,trial,status,subfail,"
          "iterations,cpu,fstar,nevals,ncorrect,variant\n")
with open(os.path.join(HERE, "results", "benchmark.csv"), "w") as f:
    f.write(header)
    for name in sorted(os.listdir(RAW)):
        label = name[:-4].split("_")[-1]
        for line in open(os.path.join(RAW, name)):
            f.write(line.rstrip("\n") + "," + label + "\n")
