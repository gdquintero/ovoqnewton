#!/bin/bash
# Gauss1_10, start 1, shifted Hessian with M = 1: the ten trials run as
# independent jobs (first=k reproduces the k-th starting point of the sequence).
cd "$(dirname "$0")"
cores=(0 2 4 6 8 10 12 14)
for k in $(seq 1 10); do
  c=${cores[$(( (k-1) % 8 ))]}
  taskset -c $c src/ovo_bench inst=$PWD/instances/Gauss1_10.txt start=1 ntrials=$k first=$k \
    seed=124456 maxit=10000 maxtime=1800 project=1 pscale=rel pert=0.25 B=H M=1 \
    out=$PWD/results/split/trial$k.csv &
  if (( k == 8 )); then wait; fi
done
wait
cat $(for k in $(seq 1 10); do echo results/split/trial$k.csv; done) > results/raw/Gauss1_10_s1_H1.csv
echo done
