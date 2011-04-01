#! /bin/sh

set -e

openmpirun -n 3 ./test_mrc_domain_multi --npx 2 --npy 2 --mrc_io_type xdmf2_parallel
#openmpirun -n 2 dtruss -f -t write ./test_mrc_domain_multi --npx 2 --npy 2 --mrc_io_type xdmf2_parallel

TEST=3
for a in reference_results/$TEST/*.xdmf; do 
    b=`basename $a`
    diff -u $a $b
done

for a in reference_results/$TEST/*.h5.dump; do 
    b=`basename $a .dump`
    h5dump $b > $b.dump
    diff -u $a $b.dump
    rm $b.dump
done

