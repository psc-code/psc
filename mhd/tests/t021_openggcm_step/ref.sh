#! /bin/bash

mpirun -n 6 ./openggcm \
    --ggcm_mhd_dipole_type float \
    --ggcm_mhd_ic_type mirdip_float \
    --xggcm_mhd_crds_type fortran \
    --ggcm_mhd_step_type fortran \
    --ggcm_mhd_satout_type none \
    --ggcm_mhd_bnd_type inoutflow_sc_ggcm_float \
    --ggcm_mhd_bnd_do_legacy \
    --ggcm_mhd_bndsw_type fortran \
    --ggcm_mhd_iono_type c_sc_ggcm_float \
    --ggcm_mhd_step_debug_dump true \
    --xmrc_io_type xdmf2 \
    --xmrc_io_sw 2 \
    \
    2>&1 | tee log
