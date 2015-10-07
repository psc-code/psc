#!/bin/sh

openmpirun -n 4 mhd_double_tearing  \
    --ggcm_mhd_ic_type double_tearing --lmx 0 --lmy 0 --lmz 0 --eta 0.0 \
    --mrc_crds_lx -2.0 --mrc_crds_ly -2.0 --mrc_crds_hx 2. --mrc_crds_hy 2. \
    --mrc_crds_lz -0.02 --mrc_crds_hz 0.02 --mrc_ts_output_every_time 0.02 \
    --mx 64 --my 64 --mz 2 --d_i 1e-4 --mrc_domain_npx 2 --mrc_domain_npy 2 \
    --mrc_ts_max_time 0.005 --mrc_ts_type rk2 --mrc_ts_dt 2e-3 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --gc_x0 0.25 --gc_x1 0.75 --gc_w 0.1 --gc_r 4.0 

exit 0
