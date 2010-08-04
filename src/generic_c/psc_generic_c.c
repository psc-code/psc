
#include "psc_generic_c.h"

struct psc_ops psc_ops_generic_c = {
  .name                   = "generic_c",
  .push_part_xz           = genc_push_part_xz,
  .push_part_yz_a         = genc_push_part_yz_a,
  .push_part_yz_b         = genc_push_part_yz_b,
};
