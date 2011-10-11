# test "cuda_1vb" against "1vb"

./test_part_yz \
    --ref_type 1vb --type cuda_1vb --moments 1st \
    --eps_particles 1e-6 --eps_fields 2e-4
