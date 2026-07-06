#!/bin/bash
rm -f salida.txt salida2.txt param.txt
touch salida.txt salida2.txt param.txt

export ALGENCAN=$HOME/algencan-3.1.1

rm -f osborne *.mod *.o

# Compilar dependencias con el gfortran de ESTA máquina
bash sort.sh   || exit 1
bash subset.sh || exit 1

# Enlazar el principal
gfortran -O3 -w -fcheck=all -g osborne.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -llapack -o osborne

delta=1.0d-4
sigmin=1.0d-1
gamma=5.d0

echo "=== Primer orden (B=0, cauchy_flag=1) ===" >> salida.txt
for ((noutliers=0; noutliers<=6; noutliers+=1))
  do
    echo $delta $sigmin $gamma $noutliers 1 > param.txt
    ./osborne >> salida.txt
  done

echo "=== Quasi-Newton (M=1, cauchy_flag=0) ===" >> salida2.txt
for ((noutliers=0; noutliers<=6; noutliers+=1))
  do
    echo $delta $sigmin $gamma $noutliers 0 > param.txt
    ./osborne >> salida2.txt
  done

echo End of execution of osborne.f90.
