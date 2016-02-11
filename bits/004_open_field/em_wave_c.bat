#! /bin/bash

# wave propagating in z direction, conducting b.c.

mpirun -n 1 ~/src/psc/src/psc_em_wave \
    --output_fields e,h,j \
    --ky 0. --kz -1. \
    --gdims_y 2 --gdims_z 32 \
    --write_tfield no \
    --write_pfield yes --pfield_step 1 \
    --stats_every 10 \
    --particles_base double \
    --fields_base c \
    --bnd_field_lo_z open \
    --bnd_field_hi_z open \
    --psc_push_particles_type 1vbec_double \
    --psc_bnd_particles_type double \
    --psc_push_fields_type c \
    --psc_bnd_type c \
    --psc_bnd_fields_type c \
    2>&1 | tee log

