
#include "psc_generic_c.h"

struct psc_ops psc_ops_generic_c = {
  .name                   = "generic_c",
  .push_part_xyz          = genc_push_part_xyz,
  .push_part_xy           = genc_push_part_xy,
  .push_part_xz           = genc_push_part_xz,
  .push_part_yz           = genc_push_part_yz,
  .push_part_z            = genc_push_part_z,
  .push_part_yz_a         = genc_push_part_yz_a,
  .push_part_yz_b         = genc_push_part_yz_b,
};

