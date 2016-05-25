#!/bin/sh 

mpirun -n 1 ../mhd_ot \
    --mrc_crds_lx 0. --mrc_crds_hx 1.0 \
    --mrc_crds_ly 0. --mrc_crds_hy 1.0 \
    --mrc_crds_lz -0.01 --mrc_crds_hz 0.01 \
    \
    --mrc_domain_mx 64 --mrc_domain_my 64 \
    --mrc_domain_npx 1 --mrc_domain_npy 1 \
    \
    --mrc_ts_output_every_time 0.01  \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    \
    --mrc_ts_max_time 1.0 \
    --ggcm_mhd_step_type c3_double \
    --xggcm_mhd_step_debug_dump \
    --ggcm_mhd_step_do_nwst \
    \
    --timelo 1000. \
    \
    2>&1 | tee log

./plot.py
