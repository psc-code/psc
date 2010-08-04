
#include "psc_generic_c.h"

struct psc_ops psc_ops_generic_c = {
  .name                   = "generic_c",
  .push_part_xz           = genc_push_part_xz,
  .push_part_yz           = genc_push_part_yz,
  .push_part_z            = genc_push_part_z,
  .push_part_yz_a         = genc_push_part_yz_a,
  .push_part_yz_b         = genc_push_part_yz_b,
};

struct psc_moment_ops psc_moment_ops_generic_c = {
  .name                   = "generic_c",
  .calc_densities         = genc_calc_densities,
};
