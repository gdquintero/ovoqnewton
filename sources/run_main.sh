rm -f salida.txt
touch salida.txt

export ALGENCAN=/opt/algencan-3.1.1
rm -f main
gfortran -O3 -w -fcheck=all -g main.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o main

# delta=0.001
# sigmin=0.1
# gamma=5

delta=1.0d-3
sigmin=1.0d-1
gamma=5.d0

echo $delta $sigmin $gamma $noutliers  > param.txt
./main >> salida.txt
