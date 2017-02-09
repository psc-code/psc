#!/bin/sh 

mpirun -n 1 ../mhd_harris \
    --mrc_domain_mx 128 --mrc_domain_my 64 --mrc_domain_mz 1 \
    --mrc_domain_npx 2 --mrc_domain_npy 4 \
    \
    --mrc_ts_output_every_time 1  \
    --ggcm_mhd_diag_fields gkeyll_e:gkeyll_i:gkeyll_em:j:rr:pp:v \
    \
    --ggcm_mhd_gk_nr_fluids 2 \
    --ggcm_mhd_gk_nr_moments 5 \
    --ggcm_mhd_gk_norm_mi_over_me 25 \
    --ggcm_mhd_gk_norm_ppi_over_ppe 5 \
    --ggcm_mhd_d_i 1. \
    \
    --mrc_ts_max_time 40.0 \
    --ggcm_mhd_step_type gkeyll \
    --ggcm_mhd_step_script ../gkeyll_step.lua \
    --xggcm_mhd_step_debug_dump \
    --ggcm_mhd_step_legacy_dt_handling false \
    \
    \
    --xtimelo 1000. \
    \
    2>&1 | tee log


./plot.py
