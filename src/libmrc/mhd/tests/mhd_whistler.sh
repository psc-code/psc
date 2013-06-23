
openmpirun -n 2 mhd_whistler \
    --mrc_crds_lx 0.0 --mrc_crds_ly 0.0 --mrc_crds_hx 0.1 --mrc_crds_hy 0.1 \
    --mrc_crds_lz 0.0 --mrc_crds_hz 2.0 \
    --init whistler --every_time 0.005 --eta 0.0 --lmx 0 --lmy 0 --lmz 0 \
    --mx 4 --my 4 --mz 128  --mrc_domain_npx 1 --mrc_domain_npz 2 \
    --mrc_ts_max_time 2.0 --mrc_ts_type rk2 --mrc_ts_dt 1e-4 \
    --pert 1e-5 --Boz 1.0 --d_i 0.25 --n0 24.0 --lambda 4 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:divb:j 
