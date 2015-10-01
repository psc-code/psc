#!/bin/sh 

mpirun -n 1 mhd_mirdip \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --lmx 0 --lmy 0 --lmz 0 \
    --mx 64 --my 64 --mz 64 \
    --mrc_ts_output_every_time 0.01  \
    --mrc_ts_max_time 1.0 \
    --do_nwst \
    --rr 6 --vx 450 --pp 6 --bz -5 \
    2>&1 | tee log

exit 0
