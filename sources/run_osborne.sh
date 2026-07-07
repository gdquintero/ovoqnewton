#!/bin/bash
rm -f salida_osborne.txt param.txt
touch salida_osborne.txt param.txt

export ALGENCAN=$HOME/algencan-3.1.1

rm -f osborne *.mod *.o

# Compilar dependencias con el gfortran de ESTA maquina
bash sort.sh   || exit 1
bash subset.sh || exit 1

# Enlazar el principal
gfortran -O3 -w -fcheck=all -g osborne.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -llapack -o osborne || exit 1

delta=1.0d-1
sigmin=1.0d-1
gamma=5.d0

# Sizes to sweep (pass as arguments, default 100)
sizes=${@:-100}

echo "method &     ndata &    obar &       o &     f(x*)  &  #it & #fcnt &      Time  &  Time/#fcnt \\\\" | tee -a salida_osborne.txt

for m in $sizes
do
  for cauchy in 1 0   # 1 = first-order (B=0), 0 = quasi-Newton (M=1)
  do
    echo "$delta $sigmin $gamma $cauchy $m" > param.txt
    ./osborne | tee -a salida_osborne.txt
    # keep a copy of the (o,f) curve tagged by size
    if [ "$cauchy" -eq 1 ]; then tag=fo; else tag=qn; fi
    cp ../output/osborne_${tag}_curve.txt ../output/osborne_${tag}_curve_${m}.txt
  done
done

echo "End of execution of osborne.f90."
