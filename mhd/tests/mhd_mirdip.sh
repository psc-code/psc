#!/bin/sh 

mpirun -n 4 mhd_mirdip \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --lmx 0 --lmy 0 --lmz 0 \
    --mx 64 --my 64 --mz 64 \
    --npx 1 --npz 4 \
    --mrc_ts_output_every_time 1.  \
    --mrc_ts_max_time 10000.0 \
    --do_nwst \
    --rr 6. --pp 6. --vx 450. --bz -5. \
    --bnd1_type sphere_sc_double \
    --bnd1_rr 3. --bnd1_pp 3. --bnd1_vx 0. --bnd1_radius 3. \
    --isphere 3. \
    --do_badval_checks false \
    --magdiffu const --diffconstant 0. \
    --xxmir 0. \
    2>&1 | tee log

exit 0
