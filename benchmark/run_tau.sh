#!/bin/bash
# Comment M3: effect of the shift margin tau in (23), Bard experiment, M = 1.
# Counts only (run on efficiency cores 16-23; CPU times are not used).
cd "$(dirname "$0")"
i=0
for tau in 0 1e-8 1e-4 1e-2; do
  for o in 0 1 2 3 4 5 6; do
    core=$((16 + i % 8)); i=$((i+1))
    taskset -c $core src/ovo_bench inst=instances/Bard_ref.txt B=H M=1 tau=$tau o=$o \
      ntrials=100 seed=123456 pscale=abs pert=0.5 project=1 \
      out=results/tau/tau${tau}_o${o}.csv &
    if (( i % 8 == 0 )); then wait; fi
  done
done
wait
echo done
