#! /bin/bash

# wave propagating in y-z plane, open b.c.

mpirun -n 1 ~/src/psc/src/psc_em_wave \
    --output_fields e,h,j \
    --k 1.41421356237309 --theta 45. \
    --amplitude_s 0. --amplitude_p 1. \
    --gdims_y 64 --gdims_z 64 \
    --nmax 201 \
    --bnd_field_lo_y open \
    --bnd_field_hi_y open \
    --bnd_field_lo_z open \
    --bnd_field_hi_z open \
    --write_tfield no \
    --write_pfield yes --pfield_step 1 \
    --stats_every 10 \
    --particles_base double \
    --fields_base c \
    --psc_push_particles_type 1vbec_double \
    --psc_bnd_particles_type double \
    --psc_push_fields_type c \
    --psc_bnd_type c \
    --psc_bnd_fields_type c \
    2>&1 | tee log

