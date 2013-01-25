# test "cuda_1vb" against "1vb"

# particle order is different, so we can't compare particles directly,
# however comparing densities works out.

./test_part_yz \
    --ref_type 1vb --type 1vb_4x4_cuda --moments 1st \
    --eps_particles 2e-6 --check_particles false \
    --eps_fields 2e-4
