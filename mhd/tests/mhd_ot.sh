openmpirun -n 2 mhd_ot \
    --mrc_crds_lx 0.0 --mrc_crds_ly 0.0 \
    --mrc_crds_hx 1.0 --mrc_crds_hy 1.0 \
    --mrc_crds_lz 0.0 --mrc_crds_hz 0.1 \
    --ggcm_mhd_diag_fields rr1:uu1:rv1:j:b1:divb:rr:pp:v:b \
    --lmx 0 --lmy 0 --lmz 0 --d_i 0.00 \
    --mx 128 --my 128 --mz 2 --npx 1 --npy 2 \
    --mrc_ts_output_every_time 0.01  \
    --mrc_ts_max_time 1.5 --mrc_ts_type rk2 --mrc_ts_dt 0.00125
