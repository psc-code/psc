#!/bin/sh 

 openmpirun -n 2 mhd_kh \
    --mrc_crds_lx -0.5 --mrc_crds_ly -0.5 \
    --mrc_crds_hx 0.5 --mrc_crds_hy 0.5 \
    --mrc_crds_lz -0.01 --mrc_crds_hz 0.01 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --lmx 0 --lmy 0 --lmz 0 \
    --mx 64 --my 64 --mz 2 --npx 1 --npy 2 \
    --mrc_ts_output_every_time 0.05 --gamma 1.4 \
    --mrc_ts_max_time 3.0 --mrc_ts_type rk2 --mrc_ts_dt 0.0025
