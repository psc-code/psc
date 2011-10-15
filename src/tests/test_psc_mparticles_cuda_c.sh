# test transfer from base "cuda" to "c" and back

./test_psc_mparticles --particles_base cuda --particles_base_flags 4096 \
    --type c --eps_particles 1e-14

