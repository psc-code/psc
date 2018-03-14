
#include "psc_heating_spot_private.h"

// ======================================================================
// psc_heating_spot

// ----------------------------------------------------------------------
// psc_heating_spot_get_H

double
psc_heating_spot_get_H(struct psc_heating_spot *spot, const double *xx)
{
  struct psc_heating_spot_ops *ops = psc_heating_spot_ops(spot);

  assert(ops && ops->get_H);
  return ops->get_H(spot, xx);
}

// ----------------------------------------------------------------------
// psc_heating_spot class

struct mrc_class_psc_heating_spot_ : mrc_class_psc_heating_spot {
  mrc_class_psc_heating_spot_() {
    name             = "psc_heating_spot";
    size             = sizeof(struct psc_heating_spot);
  }
} mrc_class_psc_heating_spot;


