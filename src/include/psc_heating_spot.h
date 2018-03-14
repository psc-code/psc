
#ifndef PSC_HEATING_SPOT_H
#define PSC_HEATING_SPOT_H

#include <psc.h>

MRC_CLASS_DECLARE(psc_heating_spot, struct psc_heating_spot);

double psc_heating_spot_get_H(struct psc_heating_spot *spot, const double *xx);

#endif

