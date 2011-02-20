
#include "psc.h"

// ======================================================================
// "none" variant of various modules, used to turn the corresponding
// functionality off.

static void
do_nothing(void)
{
}

// particle pushing

struct psc_ops psc_ops_none = {
  .name                   = "none",
  .push_part_xz           = do_nothing,
  .push_part_yz           = do_nothing,
  .push_part_z            = do_nothing,
  .push_part_yz_a         = do_nothing,
  .push_part_yz_b         = do_nothing,
};

// field advance

static void
none_push_field(struct psc_mfields *flds)
{
}

struct psc_push_field_ops psc_push_field_ops_none = {
  .name         = "none",
  .push_field_a = none_push_field,
  .push_field_b = none_push_field,
};

