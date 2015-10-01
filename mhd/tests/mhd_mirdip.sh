#!/bin/sh 

mpirun -n 1 mhd_mirdip \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --lmx 0 --lmy 0 --lmz 0 \
    --mx 64 --my 64 --mz 64 \
    --mrc_ts_output_every_time 10.  \
    --mrc_ts_max_time 1000.0 \
    --do_nwst \
    --ggcm_mhd_ic_dipole_momentx 0. \
    --ggcm_mhd_ic_dipole_momenty 0. \
    --ggcm_mhd_ic_dipole_momentz 0. \
    --rr 6 --vx 450 --pp 6 --bz 0. \
    2>&1 | tee log

exit 0
