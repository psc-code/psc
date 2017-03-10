#!/bin/sh 

RG=2634100.
dipolestrength=711.29

mpirun -n 8 ~/src/openggcm/target-build/csrc/libmrc/mhd/tests/mhd_mirdip \
    --mrc_crds_type rectilinear \
    --mrc_crds_lx -5. --mrc_crds_ly -5. --mrc_crds_lz -5. \
    --mrc_crds_hx  5. --mrc_crds_hy  5. --mrc_crds_hz  5. \
    --mrc_domain_mx 64 --mrc_domain_my 64 --mrc_domain_mz 64 \
    --mrc_domain_npx 2 --mrc_domain_npy 2 --mrc_domain_npz 2 \
    --crds_gen_x_type ggcm_yz \
    --crds_gen_x_center_spacing 0.1 \
    --crds_gen_y_type ggcm_yz \
    --crds_gen_y_center_spacing 0.1 \
    --crds_gen_z_type ggcm_yz \
    --crds_gen_z_center_spacing 0.1 \
    \
    --bnd1_type sphere_fc_cc_double \
    --bnd1_rr 56. --bnd1_pp 3800. --bnd1_radius 1. \
    --bnd1_radial_velocity 1 \
    \
    --ggcm_mhd_bnd_type inoutflow_fc_cc_double \
    --ggcm_mhd_bnd_rr 56. --ggcm_mhd_bnd_pp 3800. --ggcm_mhd_bnd_vx 140. --ggcm_mhd_bnd_bz -77. \
    \
    --ggcm_mhd_ic_type mirdip_double \
    --ggcm_mhd_ic_rr  56. --ggcm_mhd_ic_pp  3800. --ggcm_mhd_ic_vx  140. --ggcm_mhd_ic_bz  -77. \
    --ggcm_mhd_ic_rrini 56. --ggcm_mhd_ic_xmir -2.578125 --ggcm_mhd_ic_prat 1. \
    --ggcm_mhd_ic_xxx1 2.5 --ggcm_mhd_ic_xxx2 2. --ggcm_mhd_ic_stretch_tail 1. \
    --ggcm_mhd_ic_dipole_momentx 0. --ggcm_mhd_ic_dipole_momenty 0. --ggcm_mhd_ic_dipole_momentz -1. \
    --ggcm_mhd_dipole_r1lim .5 \
    \
    --ggcm_mhd_step_type mhdcc_double \
    --ggcm_mhd_step_limiter gminmod \
    --ggcm_mhd_step_riemann hll \
    --ggcm_mhd_step_bc_reconstruct \
    --ggcm_mhd_step_background \
    --ggcm_mhd_step_time_integrator euler \
    --ggcm_mhd_step_do_nwst \
    --ggcm_mhd_step_do_badval_checks false \
    --ggcm_mhd_step_divb glm \
    --ggcm_mhd_step_legacy_dt_handling no \
    --xggcm_mhd_step_debug_dump \
    --mrc_ts_output_every_time 1. \
    --mrc_ts_max_time 10000.0 \
    \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:b1:j:divb:rr:pp:v:b \
    \
    --ggcm_mhd_thx .4 \
    --ggcm_mhd_d_i 0. \
    --ggcm_mhd_magdiffu const --ggcm_mhd_diffconstant 0e4. --ggcm_mhd_diff_obnd 4 --ggcm_mhd_diffsphere .8 \
    \
    --ggcm_mhd_norm_length 1. \
    --ggcm_mhd_norm_B 1. \
    --ggcm_mhd_norm_density 1. \
    --ggcm_mhd_norm_mu0 1. \
    --ggcm_mhd_amu 1. \
    --ggcm_mhd_mu0_code 1.2566370E-06 \
    --ggcm_mhd_xxnorm0 ${RG} \
    --ggcm_mhd_rrnorm0 1.6605655E-21 \
    --ggcm_mhd_dipole_moment_distance ${RG} \
    --ggcm_mhd_dipole_dipolestrength ${dipolestrength} \
    \
    2>&1 | tee log
