#! /bin/bash

# wave propagating in y-z plane, periodic b.c.

mpirun -n 1 ~/src/psc/src/psc_em_wave \
    --output_fields e,h,j \
    --k 2.23606797749979 --theta 26.56505117707799 \
    --amplitude_s 0. --amplitude_p 1. \
    --gdims_y 32 --gdims_z 32 \
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

