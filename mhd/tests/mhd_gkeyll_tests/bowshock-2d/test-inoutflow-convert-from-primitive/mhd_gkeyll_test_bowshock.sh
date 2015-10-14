#!/bin/sh 

mpirun -n 1 ../../../mhd_gkeyll \
    --ggcm_mhd_step_debug_dump no \
    --mrc_crds_lx -1.5  --mrc_crds_hx 1.5 \
    --mrc_crds_ly -1.5  --mrc_crds_hy 1.5 \
    --mrc_crds_lz 0.  --mrc_crds_hz 1. \
    --mrc_domain_mx 128 --mrc_domain_bcx none \
    --mrc_domain_my 128 --mrc_domain_bcy none \
    --mrc_domain_mz 1   --mrc_domain_bcz periodic \
    --mrc_domain_npx 2 --mrc_domain_npy 2 --mrc_domain_npz 1 \
    \
    --ggcm_mhd_bnd_type inoutflow_gkeyll \
    --ggcm_mhd_rre 0.01 \
    --ggcm_mhd_ppe 0.0015 \
    --ggcm_mhd_vxe 1.0 \
    --ggcm_mhd_bx 0 \
    --ggcm_mhd_by 0.001 \
    --ggcm_mhd_bz 0 \
    \
    --ggcm_mhd_step_type gkeyll \
    --ggcm_mhd_step_script ../../step.lua \
    --ggcm_mhd_step_script_common ./common.lua \
    --ggcm_mhd_ic_type gkeyll \
    --ggcm_mhd_ic_script ../../init.lua \
    --ggcm_mhd_ic_script_common ./common.lua \
    \
    --mrc_ts_dt 0.1 \
    --mrc_ts_max_time 5. \
    --mrc_ts_output_every_time .25  \
    --xmrc_ts_output_every_steps 100  \
    --ggcm_mhd_diag_fields gkeyll_e:gkeyll_i:gkeyll_em \
    2>&1 | tee log

