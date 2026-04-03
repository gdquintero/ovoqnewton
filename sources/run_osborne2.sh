rm -f salida.txt param.txt
touch salida.txt param.txt

export ALGENCAN=/opt/algencan-3.1.1
rm -f osborne2
gfortran -O3 -w -fcheck=all -g osborne2.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o osborne2

# for delta in 0.5 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001
# do
#   for sigmin in 0.0001 0.001 0.01 0.1 1.0
#   do
#     for gamma in 2 5 10
#     do
#       echo 'esta=============================' >> salida.txt
#       echo 'esta=====',$delta $sigmin $gamma  >> salida.txt
#       echo 'esta=============================' >> salida.txt
#       for noutliers in {0..15}
#       do
#         echo $delta $sigmin $gamma $noutliers > param.txt
#         ./osborne2 >> salida.txt
#       done
#     done
#   done
# done

delta=1.0d-3
sigmin=1.0d-1
gamma=5.d0

for noutliers in {0..15}
  do
    echo $delta $sigmin $gamma $noutliers  > param.txt
    ./osborne2 >> salida.txt
  done


