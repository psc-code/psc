rm *.h5 *.xdmf
openmpirun -n 1 mhd_double_tearing  \
    --ggcm_mhd_ic_type double_tearing --dens0 100.0  --Bo 10.0 --pert 0.025 \
    --eps 0.3 --lambda 0.5 --lmx 0 --lmy 0 --lmz 0 --eta 0.0 \
    --mrc_crds_lx -2.0 --mrc_crds_ly -2.0 --mrc_crds_hx 2. --mrc_crds_hy 2. \
    --mrc_crds_lz -0.01 --mrc_crds_hz 0.01 --mrc_ts_output_every_time 0.01 \
    --mx 64 --my 64 --mz 2 --d_i 0.01 --mrc_domain_npx 1 --mrc_domain_npy 1 \
    --mrc_ts_max_time 8.00 --mrc_ts_type rk2 --mrc_ts_dt 1e-3 \
    --ggcm_mhd_diag_fields rr1:rv1:uu1:b1:divb # \
#    --gc_rx 2.0 --gc_x0 0.5 --gc_wx 0.25 \
#    --gc_ry 8.0 --gc_y0 0.5 --gc_wy 0.5 

paraview --state=dbg8.pvsm 
 
