# test "cuda_1st" against "1st"

./test_part_yz \
    --ref_type 1st --type cuda_1st \
    --eps_particles 1e-14 --eps_fields 1e-12
