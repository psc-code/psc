#!/bin/sh 

# test AMR functionality on an actual (made-up) two-level AMR domain
# in particular, this monitors conservation to check for proper flux correction
# at coarse/fine boundaries

DIR=`dirname $0`

openmpirun -n 1 $DIR/mhd_kh \
    --mrc_crds_lx -0.5  --mrc_crds_hx 0.5 \
    --mrc_crds_ly -0.5  --mrc_crds_hy 0.5 \
    --mrc_crds_lz -0.01 --mrc_crds_hz 0.01 \
    --pert 1e-1 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --ggcm_mhd_diag_run kh \
    --mrc_io_type xdmf2 \
    --ggcm_mhd_amr 5 \
    --mrc_domain_type amr \
    --mrc_domain_bcx periodic --mrc_domain_bcy periodic --mrc_domain_bcz none \
    --mrc_crds_type amr_uniform \
    --mx 16 --my 16 --mz 1 \
    --mrc_ts_output_every_time 0.1 --gamma 1.4 \
    --mrc_ts_max_time 2.0 \
    --mrc_ts_type step \
    --ggcm_mhd_monitor_conservation \
    --mrc_ts_conservation_every_time 0.01 \
    --ggcm_mhd_step_type vl \
    --magdiffu const \
    --do_nwst \
    2>&1 | tee log
