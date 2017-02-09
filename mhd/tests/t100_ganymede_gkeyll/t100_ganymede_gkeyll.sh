#!/bin/sh 

mpirun -n 8 ../mhd_mirdip \
    --mrc_crds_type rectilinear \
    --mrc_crds_lx -5. --mrc_crds_ly -5. --mrc_crds_lz -5. \
    --mrc_crds_hx  5. --mrc_crds_hy  5. --mrc_crds_hz  5. \
    --mrc_domain_mx 96 --mrc_domain_my 96 --mrc_domain_mz 96 \
    --mrc_domain_npx 2 --mrc_domain_npy 2 --mrc_domain_npz 2 \
    --xcrds_gen_x_type ggcm_yz \
    --crds_gen_x_center_spacing 0.08 \
    --xcrds_gen_y_type ggcm_yz \
    --crds_gen_y_center_spacing 0.08 \
    --xcrds_gen_z_type ggcm_yz \
    --crds_gen_z_center_spacing 0.08 \
    \
    --bnd1_type sphere_gkeyll \
    --bnd1_rr 56. --bnd1_pp 3800. --bnd1_radius 1. \
    --bnd1_radial_velocity 1 \
    --bnd1_test 3 \
    \
    --ggcm_mhd_bnd_type inoutflow_gkeyll \
    --ggcm_mhd_bnd_rr 56. --ggcm_mhd_bnd_pp 3800. --ggcm_mhd_bnd_vx 140. --ggcm_mhd_bnd_bz -77. \
    \
    --ggcm_mhd_ic_type mirdip_double \
    --ggcm_mhd_ic_rr  56. --ggcm_mhd_ic_pp  3800. --ggcm_mhd_ic_vx  140. --ggcm_mhd_ic_bz  -77. \
    --ggcm_mhd_ic_rrini 56. --ggcm_mhd_ic_xmir -2.578125 --ggcm_mhd_ic_prat 1. \
    --ggcm_mhd_ic_xxx1 2.5 --ggcm_mhd_ic_xxx2 2. --ggcm_mhd_ic_stretch_tail 1. \
    --ggcm_mhd_ic_dipole_momentx 0. --ggcm_mhd_ic_dipole_momenty 0. --ggcm_mhd_ic_dipole_momentz -720 \
    --ggcm_mhd_dipole_r1lim .5 \
    \
    --ggcm_mhd_step_type gkeyll \
    --ggcm_mhd_step_script ../gkeyll_step.lua \
    --ggcm_mhd_gk_nr_fluids 2 \
    --ggcm_mhd_gk_nr_moments 5 \
    --ggcm_mhd_gk_norm \
    --ggcm_mhd_gk_norm_speed_of_light 1.2 \
    --ggcm_mhd_gk_norm_mi_over_me 25. \
    --ggcm_mhd_gk_norm_ppi_over_ppe 5 \
    --ggcm_mhd_gk_norm_rr 1.4 \
    --ggcm_mhd_d_i 0.3 \
    \
    --ggcm_mhd_step_limiter gminmod \
    --ggcm_mhd_step_riemann hll \
    --ggcm_mhd_step_bc_reconstruct \
    --ggcm_mhd_step_background \
    --ggcm_mhd_step_has_ymask \
    --ggcm_mhd_step_time_integrator euler \
    --ggcm_mhd_step_do_nwst \
    --ggcm_mhd_step_do_badval_checks false \
    --ggcm_mhd_step_divb glm \
    --ggcm_mhd_step_legacy_dt_handling no \
    --xggcm_mhd_step_debug_dump \
    --mrc_ts_output_every_time 1. \
    --mrc_ts_max_time 10000.0 \
    \
    --ggcm_mhd_diag_fields gkeyll_e:gkeyll_i:gkeyll_em \
    \
    --ggcm_mhd_thx .9 \
    --ggcm_mhd_d_i 0. \
    --ggcm_mhd_magdiffu const --ggcm_mhd_diffconstant 0e4. --ggcm_mhd_diff_obnd 4 --ggcm_mhd_diffsphere .8 \
    \
    --ggcm_mhd_norm_length 2634100. \
    --ggcm_mhd_norm_B 711.29e-9 \
    --ggcm_mhd_norm_density 40e6 \
    --ggcm_mhd_earth_mag_moment 0.13725e21 \
    \
    2>&1 | tee log
