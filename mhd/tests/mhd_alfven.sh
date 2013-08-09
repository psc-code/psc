#!/bin/sh

# Alfven wave
pdim='x'

if [ $pdim = 'x' ]; then 
openmpirun -n 2 mhd_alfven \
    --mrc_crds_lx 0.0 --mrc_crds_ly 0.0 --mrc_crds_hx 1.0 --mrc_crds_hy 0.05 \
    --mrc_crds_lz 0.0 --mrc_crds_hz 0.05 \
    --init bw  --ggcm_mhd_ic_pdim x --lmx 0 --lmy 0 --lmz 0 \
    --mx 1024 --my 2 --mz 2 --mrc_ts_dt 0.001 --npx 2 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:j \
    --mrc_ts_max_time 2.0 --mrc_ts_output_every_time 0.02
fi 


if [ $pdim = 'xx' ]; then 
openmpirun -n 2 mhd_bw \
    --mrc_crds_lx 0.0 --mrc_crds_ly 0.0 --mrc_crds_hx 1.0 --mrc_crds_hy 0.05 \
    --mrc_crds_lz 0.0 --mrc_crds_hz 0.05 --gamma 2.0 --bcx NONE \
    --init bw  --ggcm_mhd_ic_pdim x --lmx 0 --lmy 0 --lmz 0 \
    --mx 512 --my 2 --mz 2 --mrc_ts_dt 0.00025 --npx 2 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:j \
    --mrc_ts_max_time 0.5 --mrc_ts_output_every_time 0.002
fi 



if [ $pdim = 'xx256' ]; then 
openmpirun -n 1 mhd_bw \
    --mrc_crds_lx 0.0 --mrc_crds_ly 0.0 --mrc_crds_hx 1.0 --mrc_crds_hy 0.2 \
    --mrc_crds_lz 0.0 --mrc_crds_hz 0.2 --gamma 2.0 --bcx NONE \
    --init bw  --ggcm_mhd_ic_pdim x --lmx 0 --lmy 0 --lmz 0 \
    --mx 256 --my 4 --mz 4 --mrc_ts_dt 0.00025 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:divb:j \
    --mrc_ts_max_time 0.5 --mrc_ts_output_every_time 0.002
fi 




if [ $pdim = 'y' ]; then 
openmpirun -n 1 mhd_bw \
    --mrc_crds_lx 0.0 --mrc_crds_ly 0.0 --mrc_crds_hx 0.05 --mrc_crds_hy 1.0 \
    --mrc_crds_lz 0.0 --mrc_crds_hz 0.05 --gamma 2.0 --bcy NONE \
    --init bw  --ggcm_mhd_ic_pdim y --lmx 0 --lmy 0 --lmz 0 \
    --mx 4 --my 1024 --mz 4 --mrc_ts_dt 0.00000625 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:divb:j \
    --mrc_ts_max_time 0.25 --mrc_ts_output_every_time 0.002
fi 


if [ $pdim = 'z' ]; then 
openmpirun -n 1 mhd_bw \
    --mrc_crds_lx 0. --mrc_crds_ly 0.0 --mrc_crds_hx 0.05 --mrc_crds_hy 0.05 \
    --mrc_crds_lz 0. --mrc_crds_hz 1.0 --gamma 2.0 --bcz NONE \
    --init bw  --ggcm_mhd_ic_pdim z --lmx 0 --lmy 0 --lmz 0 \
    --mx 4 --my 4 --mz 1024 --mrc_ts_dt 0.00000625 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:divb:j \
    --mrc_ts_max_time 0.25 --mrc_ts_output_every_time 0.002
fi 


exit 0 
