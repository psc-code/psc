#!/bin/sh 

mpirun -n 1  ~/src/openggcm/target-build/csrc/libmrc/mhd/tests/mhd_gkeyll \
    --mrc_crds_lx 0.  --mrc_crds_hx 6.28318530718 \
    --mrc_crds_ly 0.  --mrc_crds_hy 6.28318530718 \
    --mrc_crds_lz 0.  --mrc_crds_hz 1. \
    --mrc_domain_mx 128 --mrc_domain_bcx periodic \
    --mrc_domain_my 128 --mrc_domain_bcy periodic \
    --mrc_domain_mz 1   --mrc_domain_bcz periodic \
    --mrc_domain_npx 1 --mrc_domain_npy 1 --mrc_domain_npz 1 \
    --ggcm_mhd_bnd_type none \
    --ggcm_mhd_bndsw_type none \
    \
    --ggcm_mhd_step_type gkeyll \
    --ggcm_mhd_step_script step.lua \
    --ggcm_mhd_step_script_common common.lua \
    --ggcm_mhd_ic_type gkeyll \
    --ggcm_mhd_ic_script init.lua \
    \
    --mrc_ts_dt 0.1 \
    --mrc_ts_max_time 2.5 \
    --mrc_ts_output_every_time .1  \
    --xmrc_ts_output_every_steps 100  \
    --ggcm_mhd_diag_fields gkeyll_e:gkeyll_i:gkeyll_em \
    2>&1 | tee log

