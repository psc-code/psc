#! /bin/bash

# sample harris run
# underresolved, but we only care about continuity / gauss's law here

mpirun -n 8 ~/src/psc/src/psc_harris \
    --mi_over_me 25. \
    --gdims_y 200 --gdims_z 800 \
    --npy 8 --npz 16 \
    --nmax 51 \
    --output_fields e,h,j,n_1st_single,v_1st_single \
    --write_tfield no \
    --write_pfield yes --pfield_step 100 \
    --stats_every 10 \
    --psc_checks_continuity_every_step 10 \
    --psc_checks_continuity_verbose \
    --psc_checks_continuity_threshold 1e-5 \
    --psc_checks_continuity_dump_always \
    --psc_checks_gauss_every_step 10 \
    --psc_checks_gauss_verbose \
    --psc_checks_gauss_threshold 1e-5 \
    --psc_checks_gauss_dump_always \
    --particles_base double \
    --fields_base c \
    --psc_push_particles_type 1vbec_double \
    --psc_bnd_particles_type double \
    --psc_push_fields_type c \
    --psc_bnd_type c \
    --psc_bnd_fields_type c \
    2>&1 | tee log

