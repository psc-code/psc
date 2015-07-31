#!/bin/sh 

mpirun -n 1 ../mhd_gkeyll \
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
    --ggcm_mhd_bndsw_type constant_5m \
    --ggcm_mhd_bndsw_rre 0.00038461538461538 \
    --ggcm_mhd_bndsw_ppe 0.00075 \
    --ggcm_mhd_bndsw_vxe 1.0 \
    --ggcm_mhd_bndsw_rri 0.0096153846153846 \
    --ggcm_mhd_bndsw_ppi 0.00075 \
    --ggcm_mhd_bndsw_vxi  1.0 \
    --ggcm_mhd_bndsw_ex 0 \
    --ggcm_mhd_bndsw_ey 0 \
    --ggcm_mhd_bndsw_ez 0 \
    --ggcm_mhd_bndsw_bx 0 \
    --ggcm_mhd_bndsw_by 0.001 \
    --ggcm_mhd_bndsw_bz 0 \
    \
    --ggcm_mhd_step_type gkeyll \
    --ggcm_mhd_step_script step.lua \
    --ggcm_mhd_step_script_common common.lua \
    --ggcm_mhd_ic_type gkeyll \
    --ggcm_mhd_ic_script init.lua \
    \
    --mrc_ts_dt 0.1 \
    --mrc_ts_max_time 5. \
    --mrc_ts_output_every_time .25  \
    --xmrc_ts_output_every_steps 100  \
    --ggcm_mhd_diag_fields gkeyll_e:gkeyll_i:gkeyll_em \
    2>&1 | tee log

