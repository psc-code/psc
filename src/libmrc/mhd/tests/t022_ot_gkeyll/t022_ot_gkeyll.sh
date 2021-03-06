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
    --ggcm_mhd_diag_fields gkeyll_e:gkeyll_i:gkeyll_em:j:rr:pp:v \
    \
    --ggcm_mhd_gk_nr_fluids 2 \
    --ggcm_mhd_gk_nr_moments 5 \
    --ggcm_mhd_d_i .05 \
    \
    --mrc_ts_max_time 1.0 \
    --ggcm_mhd_step_type gkeyll \
    --ggcm_mhd_step_script ../gkeyll_step.lua \
    --xggcm_mhd_step_debug_dump \
    --ggcm_mhd_step_legacy_dt_handling false \
    \
    \
    --xtimelo 1000. \
    \
    2>&1 | tee log

#./plot.py
