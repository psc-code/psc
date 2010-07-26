
#ifndef PSC_PULSE_H
#define PSC_PULSE_H

#include "psc.h"

#include <stdlib.h>
#include <string.h>

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

static inline struct psc_pulse *
psc_pulse_create(size_t size, struct psc_pulse_ops *ops)
{
  struct psc_pulse *pulse = malloc(size);
  memset(pulse, 0, size);
  pulse->ops = ops;

  return pulse;
}


#endif
