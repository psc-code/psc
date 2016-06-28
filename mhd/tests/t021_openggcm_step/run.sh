#! /bin/bash

mpirun -n 6 ./openggcm \
    --ggcm_mhd_dipole_type float \
    --ggcm_mhd_ic_type mirdip_float \
    --xggcm_mhd_crds_type fortran \
    --ggcm_mhd_step_type c_float \
    --ggcm_mhd_satout_type none \
    --ggcm_mhd_bnd_type inoutflow_sc_ggcm_float \
    --ggcm_mhd_bnd_do_legacy \
    --ggcm_mhd_bndsw_type fortran \
    --ggcm_mhd_iono_type c_sc_ggcm_float \
    --ggcm_mhd_step_debug_dump true \
    --xmrc_io_type xdmf2 \
    --xmrc_io_sw 2 \
    \
    --mhd_primvar c \
    --mhd_primbb c \
    --mhd_zmaskn c \
    --mhd_rmaskn c \
    --mhd_newstep c \
    --mhd_push c \
    --mhd_pushpred c \
    --mhd_pushcorr c \
    --mhd_pushfluid1 c \
    --mhd_pushfield1 c \
    --mhd_pushfluid2 c \
    --mhd_pushfield2 c \
    --mhd_push_ej c \
    --mhd_pfie3 c \
    --mhd_calce c \
    --mhd_bpush1 c \
    --mhd_calc_resis c \
    \
    2>&1 | tee log
