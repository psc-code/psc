
#ifndef PSC_PULSE_H
#define PSC_PULSE_H

#include "psc.h"

struct psc_pulse;

struct psc_pulse_ops {
  const char *name;
  void (*destroy)(struct psc_pulse *);
  double (*field)(struct psc_pulse *,
		  double x, double y, double z, double t);
};  

struct psc_pulse {
  struct psc_pulse_ops *ops;
};

struct psc_pulse *psc_pulse_p_z1_short_create(void);

#endif
