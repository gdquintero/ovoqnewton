#!/bin/bash
rm -f salida_scalbard.txt param.txt
touch salida_scalbard.txt param.txt

export ALGENCAN=$HOME/algencan-3.1.1

rm -f scalbard *.mod *.o

# Compilar dependencias con el gfortran de ESTA maquina
bash sort.sh   || exit 1
bash subset.sh || exit 1

# Enlazar el principal
gfortran -O3 -w -fcheck=all -g scalbard.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -llapack -o scalbard || exit 1

delta=1.0d-1
sigmin=1.0d-1
gamma=5.d0

# Sizes to sweep (pass as arguments, default: 100 1000 10000)
sizes=${@:-"100 1000 10000"}

echo "method &     ndata &    obar &       o &     f(x*)  &  #it & #fcnt &      Time  &  Time/#fcnt \\\\" | tee -a salida_scalbard.txt

for m in $sizes
do
  # 1) Generate the synthetic data for this m (Python), 2) box-constrained LS start
  python data_scalbard.py "$m"   || exit 1
  python least_squares.py scalbard || exit 1
  # keep a per-m copy of the data and LS start (sol_ls_scalbard.txt is overwritten each m)
  cp ../data/scalbard.txt            ../data/scalbard_${m}.txt
  cp ../output/sol_ls_scalbard.txt   ../output/sol_ls_scalbard_${m}.txt

  for cauchy in 1 0   # 1 = first-order (B=0), 0 = quasi-Newton (M=1)
  do
    echo "$delta $sigmin $gamma $cauchy" > param.txt
    ./scalbard | tee -a salida_scalbard.txt
    if [ "$cauchy" -eq 1 ]; then tag=fo; else tag=qn; fi
    cp ../output/scalbard_${tag}_curve.txt ../output/scalbard_${tag}_curve_${m}.txt
  done
done

echo "End of execution of scalbard.f90."
