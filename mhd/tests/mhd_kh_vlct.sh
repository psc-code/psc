#!/bin/sh 

openmpirun -n 4 ~/src/openggcm/target-build/csrc/libmrc/mhd/tests/mhd_kh \
    --mrc_crds_lx -0.5 --mrc_crds_hx 0.5 \
    --mrc_crds_ly -0.5 --mrc_crds_hy 0.5 \
    --mrc_crds_lz -0.01 --mrc_crds_hz 0.01 \
    --pert 0e-1 --pert_random 1e-2 \
    --B0x 0 --B0z_harris .5 --lambda .02 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --ggcm_mhd_diag_run kh_vlct \
    --mx 256 --my 256 --mz 1 --npx 2 --npy 2 \
    --mrc_ts_output_every_time 0.1 --gamma 1.666667 \
    --mrc_ts_max_time 5.0 \
    --ggcm_mhd_step_type c3_double \
    --magdiffu const \
    --do_nwst \
    2>&1 | tee log
