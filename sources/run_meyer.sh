#!/bin/bash
rm -f salida.txt param.txt
touch salida.txt param.txt

export ALGENCAN=$HOME/algencan-3.1.1

rm -f meyer *.mod *.o

# Compilar dependencias con el gfortran de ESTA máquina
bash sort.sh   || exit 1
bash subset.sh || exit 1

# Enlazar el principal
gfortran -O3 -w -fcheck=all -g meyer.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -llapack -o meyer

delta=1.0d-2
sigmin=1.0d-2
gamma=2.d0

for ((noutliers=0; noutliers<=0; noutliers+=1))
  do
    echo $delta $sigmin $gamma $noutliers  > param.txt
    ./meyer >> salida.txt
  done

# for delta in 0.1 0.01 0.001 0.0001
# do
#   for sigmin in 0.01 0.1 1.0
#   do
#     for gamma in 2 5 10
#     do
#       echo 'esta=============================' >> salida.txt
#       echo 'esta=====',$delta $sigmin $gamma  >> salida.txt
#       echo 'esta=============================' >> salida.txt
#       for noutliers in {0..0}
#       do
#         echo $delta $sigmin $gamma $noutliers > param.txt
#         ./meyer >> salida.txt
#       done
#     done
#   done
# done

echo End of execution of meyer.f90.