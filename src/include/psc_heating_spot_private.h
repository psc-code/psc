
#ifndef PSC_HEATING_SPOT_PRIVATE_H
#define PSC_HEATING_SPOT_PRIVATE_H

#include <psc_heating_spot.h>

struct psc_heating_spot {
  struct mrc_obj obj;
};

struct psc_heating_spot_ops {
  MRC_SUBCLASS_OPS(struct psc_heating_spot);
  double (*get_H)(struct psc_heating_spot *spot, const double x[3]);
};

#define psc_heating_spot_ops(spot) ((struct psc_heating_spot_ops *)((spot)->obj.ops))

#endif
