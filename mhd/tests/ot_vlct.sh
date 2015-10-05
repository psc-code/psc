#!/bin/bash 

openmpirun -n 2 mhd_ot \
    --mx 128 --my 128 --mz 4 --npx 2 --npy 1 \
    --magdiffu const \
    --mrc_ts_output_every_time 0.01  \
    --mrc_ts_max_time .1 \
    --mrc_ts_dt 0.001 \
    --mrc_ts_type step \
    --ggcm_mhd_step_type vlct \
    --ggcm_mhd_diag_fields rr:pp:v:b:rr1:rv1:uu1:b1:divb \
    --ggcm_mhd_diag_run ot \
    --mhd_riemann_type hlld \
    --resnorm 53.5848e6 \
    --diffconstant 1e4 \
    --diffsphere 0. \
    --do_nwst 1 
