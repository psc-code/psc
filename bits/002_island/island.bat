#! /bin/bash
#PBS -l nodes=32:ppn=32
#PBS -j oe

cd $PBS_O_WORKDIR

# sample island coalescence run

aprun -n 1024 ~/src/psc/src/psc_island_coalescence \
    --np_y 128 --np_z 256 \
    --nicell 300 \
    --gdims_y 1024 --gdims_z 2048 \
    --nmax 100001 \
    --lambda 10. \
    --write_tfield no \
    --write_pfield yes --pfield_step 1000 \
    --output_fields e,h,j,n_1st_single,v_1st_single,divb \
    --stats_every 10 \
    --psc_balance_every 5000 \
    --fields_base single \
    --psc_push_fields_type single \
    --psc_push_fields_variant 1 \
    --psc_bnd_type single \
    --psc_bnd_fields_type single \
    --particles_base single \
    --psc_push_particles_type 1vbec_single \
    --psc_bnd_particles_type single \
    2>&1 | tee log

