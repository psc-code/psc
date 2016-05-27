#!/bin/sh 

mpirun -n 4 ../mhd_mirdip \
    --mrc_crds_type rectilinear \
    --mrc_crds_lx -3.5 --mrc_crds_ly -3.5 --mrc_crds_lz -3.5 \
    --mrc_crds_hx  3.5 --mrc_crds_hy  3.5 --mrc_crds_hz  3.5 \
    --crds_gen_x_type ggcm_yz \
    --crds_gen_x_center_spacing 0.04 \
    --crds_gen_y_type ggcm_yz \
    --crds_gen_y_center_spacing 0.04 \
    --crds_gen_z_type ggcm_yz \
    --crds_gen_z_center_spacing 0.04 \
    \
    --mrc_domain_mx 64 --mrc_domain_my 64 --mrc_domain_mz 64 \
    --mrc_domain_npx 2 --mrc_domain_npy 1 --mrc_domain_npz 2 \
    \
    --ggcm_mhd_ic_rr 56. --ggcm_mhd_ic_pp 3800. --ggcm_mhd_ic_vx 140. --ggcm_mhd_bz -77. \
    --ggcm_mhd_ic_rrini 50. \
    --ggcm_mhd_ic_prat 1. \
    --ggcm_mhd_ic_xxx1 2.5 --ggcm_mhd_ic_xxx2 2.2 \
    --ggcm_mhd_ic_xmir -2.6 \
    --ggcm_mhd_dipole_r1lim .5 \
    \
    --ggcm_mhd_bnd_type inoutflow_fc_cc_double \
    \
    --bnd1_type sphere_fc_cc_double \
    --bnd1_rr 55. --bnd1_pp 3800. --bnd1_vx 0. --bnd1_radius 1. \
    \
    --mrc_ts_output_every_time 1.  \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:b1:j:divb:rr:pp:v:b:ymask \
    \
    --mrc_ts_max_time 0. \
    --ggcm_mhd_step_type mhdcc_double \
    --ggcm_mhd_step_do_nwst \
    \
    --ggcm_mhd_norm_length 2634100. \
    --ggcm_mhd_norm_B 711.29e-9 \
    --ggcm_mhd_norm_density 40e6 \
    --ggcm_mhd_earth_mag_moment 0.13725e21 \
    \
    2>&1 | tee log

./plot.py
