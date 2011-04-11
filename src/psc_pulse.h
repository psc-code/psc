
#ifndef PSC_PULSE_H
#define PSC_PULSE_H

#include "psc.h"
#include <mrc_params.h>

#include <stdlib.h>
#include <string.h>

// ----------------------------------------------------------------------
// psc_pulse

MRC_CLASS_DECLARE(psc_pulse, struct psc_pulse);

double psc_pulse_field_s(struct psc_pulse *pulse, double x, double y, double z, double t);
double psc_pulse_field_p(struct psc_pulse *pulse, double x, double y, double z, double t);

#endif
