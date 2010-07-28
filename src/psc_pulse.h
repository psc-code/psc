
#ifndef PSC_PULSE_H
#define PSC_PULSE_H

#include "psc.h"
#include "util/params.h"

#include <stdlib.h>
#include <string.h>

// ----------------------------------------------------------------------
// psc_pulse

struct psc_pulse;

struct psc_pulse_ops {
  const char *name;
  size_t ctx_size;
  struct param *ctx_descr;
  void (*destroy)(struct psc_pulse *);
  void (*setup)(struct psc_pulse *);
  double (*field)(struct psc_pulse *,
		  double x, double y, double z, double t);
};  

struct psc_pulse {
  struct psc_pulse_ops *ops;
  bool is_setup;
  void *ctx;
};

struct psc_pulse *psc_pulse_create(struct psc_pulse_ops *ops, void *prm);
void psc_pulse_destroy(struct psc_pulse *pulse);

static inline void
psc_pulse_setup(struct psc_pulse *pulse)
{
  if (pulse->ops->setup) {
    pulse->ops->setup(pulse);
  }
  pulse->is_setup = true;
}

static inline double
psc_pulse_field(struct psc_pulse *pulse, double x, double y, double z, double t)
{
  if (!pulse->is_setup) {
    psc_pulse_setup(pulse);
  }
  return pulse->ops->field(pulse, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z1 // FIXME rename ->gauss

struct psc_p_pulse_z1_param {
  double xm, ym, zm; // location of pulse center at time 0 in m 
  double dxm, dym, dzm; // width of pulse in m
};

struct psc_pulse *psc_pulse_p_z1_short_create(struct psc_p_pulse_z1_param *prm);  // FIXME, rename -> gauss

// ----------------------------------------------------------------------
// psc_p_pulse_z1_flattop

struct psc_p_pulse_z1_flattop_param {
  double xm, ym, zm; // location of pulse center at time 0 in m 
  double dxm, dym, dzm; // slope of pulse in m
  double zb; // width of pulse in m
};

struct psc_pulse *psc_pulse_p_z1_flattop_create(struct psc_p_pulse_z1_flattop_param *prm);

#endif
