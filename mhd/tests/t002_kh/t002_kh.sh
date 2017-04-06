#!/bin/sh 

mpirun -n 1 ../mhd_kh \
    --mrc_crds_lx -0.5 --mrc_crds_ly -0.5 \
    --mrc_crds_hx 0.5 --mrc_crds_hy 0.5 \
    \
    --mrc_domain_mx 64 --mrc_domain_my 64 \
    --mrc_domain_npx 1 --mrc_domain_npy 1 \
    \
    --ggcm_mhd_ic_pert 1e-1 \
    \
    --ggcm_mhd_gamma 1.4 \
    \
    --mrc_ts_output_every_time 0.1 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    \
    --mrc_ts_max_time 3. \
    --ggcm_mhd_step_type c3_double \
    --ggcm_mhd_step_do_nwst \
    \
    2>&1 | tee log

./plot.py
