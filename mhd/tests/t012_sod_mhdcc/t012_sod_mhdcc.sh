#!/bin/sh 

mpirun -n 1 ../mhd_sod \
    \
    --mrc_domain_mx 128 \
    --mrc_domain_npx 1 \
    \
    --mrc_ts_output_every_time 0.1  \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    \
    --mrc_ts_max_time .3 \
    --ggcm_mhd_step_type mhdcc_double \
    --ggcm_mhd_step_limiter gminmod \
    --ggcm_mhd_step_riemann hll \
    --xggcm_mhd_step_debug_dump \
    --ggcm_mhd_step_do_nwst \
    \
    --xtimelo 1000. \
    \
    2>&1 | tee log

./plot.py
