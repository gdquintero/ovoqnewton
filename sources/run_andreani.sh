#!/bin/bash
rm -f salida.txt param.txt
touch salida.txt param.txt

export ALGENCAN=$HOME/algencan-3.1.1

rm -f andreani *.mod *.o

# Compilar dependencias con el gfortran de ESTA máquina
bash sort.sh   || exit 1
bash subset.sh || exit 1

# Enlazar el principal
gfortran -O3 -w -fcheck=all -g andreani.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -llapack -o andreani

delta=1.0d-1
sigmin=1.0d-1
gamma=5.d0

for ((noutliers=0; noutliers<=0; noutliers+=1))
  do
    echo $delta $sigmin $gamma $noutliers  > param.txt
    ./andreani >> salida.txt
  done

echo End of execution of andreani.f90.