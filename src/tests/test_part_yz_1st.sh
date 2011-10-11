# test "1st" against "generic_c"
# whereas the latter is 2nd order, the fields are linear,
# so the particles should see the same fields.

./test_part_yz \
    --ref_type 1st --type generic_c \
    --eps_particles 1e-14 --eps_fields 1e-12
