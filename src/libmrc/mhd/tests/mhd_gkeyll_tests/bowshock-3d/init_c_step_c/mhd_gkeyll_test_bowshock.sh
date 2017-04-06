#!/bin/sh 

mpirun -n 1 ../../../mhd_gkeyll \
    --ggcm_mhd_step_debug_dump no \
    --mrc_crds_lx -0.75  --mrc_crds_hx 0.75 \
    --mrc_crds_ly -0.75  --mrc_crds_hy 0.75 \
    --mrc_crds_lz -0.375  --mrc_crds_hz 0.375 \
    --mrc_domain_mx 64 --mrc_domain_bcx none \
    --mrc_domain_my 64 --mrc_domain_bcy none \
    --mrc_domain_mz 32 --mrc_domain_bcz none \
    --mrc_domain_npx 2 --mrc_domain_npy 2 --mrc_domain_npz 1 \
    \
    --ggcm_mhd_bnd_type inoutflow_fc_gkeyll \
    --ggcm_mhd_bnd_rr 0.01 \
    --ggcm_mhd_bnd_pp 0.0015 \
    --ggcm_mhd_bnd_vx 1.0 \
    --ggcm_mhd_bnd_ex 0 \
    --ggcm_mhd_bnd_ey 0 \
    --ggcm_mhd_bnd_ez 0 \
    --ggcm_mhd_bnd_bx 0 \
    --ggcm_mhd_bnd_by 0.001 \
    --ggcm_mhd_bnd_bz 0 \
    \
    --ggcm_mhd_step_type vlct \
    \
    --ggcm_mhd_ic_type bowshock3d \
    --ggcm_mhd_ic_x_obstacle -0.375 \
    --ggcm_mhd_ic_rho_obstacle 100 \
    --ggcm_mhd_ic_rho0 0.01 \
    --ggcm_mhd_ic_By0 0.001 \
    \
    --mrc_ts_dt 0.1 \
    --mrc_ts_max_time 5. \
    --mrc_ts_output_every_time .25  \
    --xmrc_ts_output_every_steps 100  \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:b1:e_ec \
    2>&1 | tee log

