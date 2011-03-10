
#include "psc.h"

// ======================================================================
// "none" variant of various modules, used to turn the corresponding
// functionality off.

static void
none_push_part(mfields_base_t *flds, mparticles_base_t *particles)
{
}

// particle pushing

struct psc_ops psc_ops_none = {
  .name                   = "none",
  .push_part_xz           = none_push_part,
  .push_part_yz           = none_push_part,
  .push_part_z            = none_push_part,
  .push_part_yz_a         = none_push_part,
  .push_part_yz_b         = none_push_part,
};

