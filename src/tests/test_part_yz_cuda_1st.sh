# test "cuda_1st" against "1st"

# particle order is different, so we can't compare particles directly,
# however comparing densities works out.

./test_part_yz \
    --ref_type 1st --type cuda_1st --moments 1st \
    --eps_particles 2e-6 --check_particles false \
    --eps_fields 2e-4
