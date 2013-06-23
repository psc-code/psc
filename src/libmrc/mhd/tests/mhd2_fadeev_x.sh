
openmpirun -n 1 ./mhd2  \
    --ggcm_mhd_ic_type fadeev --dens0 100.0  --Bo 10.0 --pert 0.05 \
    --eps 0.3 --lambda 0.5 --lmx 0 --lmy 0 --lmz 0 \
    --mrc_crds_lx -1. --mrc_crds_ly -4.0 --mrc_crds_hx 1. --mrc_crds_hy 4.0 \
    --mrc_crds_lz -0.01 --mrc_crds_hz 0.01 --mrc_ts_output_every_time 0.01 --eta 1e-4 \
    --mx 128 --my 128 --mz 2 --d_i 0. --mrc_domain_npx 1 --mrc_domain_npy 1 \
    --mrc_ts_max_time 6.0 --mrc_ts_type rk2 --mrc_ts_dt 5e-3 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:divb \
    --gc_rx 4.0 --gc_x0 0.5 --gc_wx 0.25 \
    --gc_ry 8.0 --gc_y0 0.5 --gc_wy 0.5 
