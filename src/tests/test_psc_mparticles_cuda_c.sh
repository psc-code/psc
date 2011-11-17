# test transfer from base "cuda" to "c" and back

./test_psc_mparticles --particles_base cuda --type c --eps_particles 1e-14 \
    --psc_push_particles_type cuda_1vb

