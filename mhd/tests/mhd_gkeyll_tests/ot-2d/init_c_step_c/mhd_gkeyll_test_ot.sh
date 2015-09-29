#!/bin/sh 

mpirun -n 1  ../../../mhd_gkeyll \
    --mrc_crds_lx 0.  --mrc_crds_hx 1. \
    --mrc_crds_ly 0.  --mrc_crds_hy 1. \
    --mrc_crds_lz 0.  --mrc_crds_hz 1. \
    --mrc_domain_mx 128 --mrc_domain_bcx periodic \
    --mrc_domain_my 128 --mrc_domain_bcy periodic \
    --mrc_domain_mz 1   --mrc_domain_bcz periodic \
    --mrc_domain_npx 2 --mrc_domain_npy 2 --mrc_domain_npz 1 \
    --ggcm_mhd_bnd_type none \
    --ggcm_mhd_bndsw_type none \
    \
    --ggcm_mhd_step_type vlct \
    --ggcm_mhd_ic_type ot \
    \
    --mrc_ts_dt 0.001 \
    --mrc_ts_max_time .51 \
    --mrc_ts_output_every_time .05  \
    --xmrc_ts_output_every_steps 100  \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:e_ec \
    2>&1 | tee log

