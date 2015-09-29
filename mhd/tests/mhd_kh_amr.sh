#!/bin/sh 

openmpirun -n 1 ./mhd_kh \
    --mrc_crds_lx -0.5 --mrc_crds_ly -0.5 \
    --mrc_crds_hx 0.5 --mrc_crds_hy 0.5 \
    --mrc_crds_lz -0.01 --mrc_crds_hz 0.01 \
    --pert 0.01 --pert_random .0 --B0 .129 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:b1:divb:rr:pp:v:b \
    --ggcm_mhd_diag_run kh \
    --mrc_io_type xdmf2 --mrc_io_sw 0 \
    --ggcm_mhd_amr 7 \
    --mrc_domain_type amr \
    --mrc_domain_bcx periodic --mrc_domain_bcy periodic --mrc_domain_bcz none \
    --mrc_crds_type amr_uniform \
    --mx 16 --my 16 --mz 1 \
    --mrc_ts_output_every_time 0.1 \
    --mrc_ts_conservation_every_time 0.01 \
    --mrc_ts_max_time 5.0 \
    --mrc_ts_type step \
    --ggcm_mhd_step_type vlct \
    --magdiffu const \
    --do_nwst \
    2>&1 | tee log
