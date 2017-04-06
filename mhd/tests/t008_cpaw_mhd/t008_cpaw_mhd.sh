#!/bin/sh 

mpirun -n 1 ../mhd_cpaw \
    --mrc_domain_mx 32 --mrc_domain_my 1 --mrc_domain_mz 1 \
    --mrc_domain_npx 1 --mrc_domain_npy 1 --mrc_domain_npz 1 \
    \
    --ggcm_mhd_d_i 0. \
    \
    --mrc_ts_max_time 2. \
    --ggcm_mhd_step_type mhdcc_double \
    --ggcm_mhd_step_do_nwst \
    --ggcm_mhd_thx .8 \
    \
    --mrc_ts_output_every_time .1 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:b1:j:divb:rr:pp:v:b \
    \
    2>&1 | tee log

./plot.py
