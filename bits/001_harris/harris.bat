#! /bin/bash

#PBS -l nodes=4:ppn=32
#PBS -j oe

cd $PBS_O_WORKDIR

# sample harris run

aprun -n 128 ~/src/psc/src/psc_harris \
    --mi_over_me 25. \
    --npy 8 --npz 16 \
    --nmax 50001 \
    --output_fields e,h,j,n_1st_single,v_1st_single \
    --write_tfield no \
    --write_pfield yes --pfield_step 100 \
    --stats_every 10 \
    --gdims_in_terms_of_cells \
    --particles_base single \
    --fields_base single \
    --psc_push_particles_type 1vbec_single \
    --psc_bnd_particles_type single \
    --psc_push_fields_type single \
    --psc_bnd_type single \
    --psc_bnd_fields_type single \
    2>&1 | tee log

