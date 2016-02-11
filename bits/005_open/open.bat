#! /bin/bash

#

mpirun -n 8 ~/src/psc/src/psc_test_open \
    --mi_over_me 25. \
    --gdims_y 64 --gdims_z 128 \
    --np_y 4 --np_z 2 \
    --output_fields e,h,j,n_1st_single,v_1st_single \
    --write_tfield no \
    --write_pfield yes --pfield_step 10 \
    --stats_every 10 \
    --particles_base single \
    --fields_base single \
    --psc_push_particles_type 1vbec_single \
    --psc_bnd_particles_type single \
    --psc_push_fields_type single \
    --psc_bnd_type single \
    --psc_bnd_fields_type single \
    2>&1 | tee log

