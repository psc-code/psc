# test "cuda_1vb" against "1vb"

# particle order is different, so we can't compare particles directly,
# however comparing densities works out.

./test_part_yz \
    --ref_type 1vb --type cuda_1vb --moments 1st \
    --eps_particles 1e-6 --check_particles false \
    --eps_fields 2e-4
