
#ifndef PSC_PULSE_H
#define PSC_PULSE_H

#include "psc.h"
#include <mrc_params.h>

#include <stdlib.h>
#include <string.h>

// ----------------------------------------------------------------------
// psc_pulse

MRC_CLASS_DECLARE(psc_pulse, struct psc_pulse);

struct psc_pulse_ops {
  MRC_SUBCLASS_OPS(struct psc_pulse);
  double (*field_s)(struct psc_pulse *,
		    double x, double y, double z, double t);
  double (*field_p)(struct psc_pulse *,
		    double x, double y, double z, double t);
};  

struct psc_pulse {
  struct mrc_obj obj;
  bool is_setup;
};

extern struct psc_pulse_ops psc_pulse_gauss_ops;
extern struct psc_pulse_ops psc_pulse_flattop_ops;

#define psc_pulse_ops(pulse) ((struct psc_pulse_ops *)((pulse)->obj.ops))

void psc_pulse_ini(struct psc_pulse *pulse, struct psc_pulse_ops *ops, void *prm);

static inline double
psc_pulse_field_s(struct psc_pulse *pulse, double x, double y, double z, double t)
{
  if (!pulse->is_setup) {
    psc_pulse_setup(pulse);
  }
  return psc_pulse_ops(pulse)->field_s(pulse, x, y, z, t);
}

static inline double
psc_pulse_field_p(struct psc_pulse *pulse, double x, double y, double z, double t)
{
  if (!pulse->is_setup) {
    psc_pulse_setup(pulse);
  }
  return psc_pulse_ops(pulse)->field_p(pulse, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_pulse_gauss

struct psc_pulse_gauss {
  double xm, ym, zm; // location of pulse center at time 0 in m 
  double dxm, dym, dzm; // width of pulse in m
  double amplitude_p;   // max amplitude, p-polarization
  double amplitude_s;   // max amplitude, s-polarization
  double phase_p;       // CEP-phase  (from -pi to pi)
  double phase_s;       // CEP-phase  (from -pi to pi)
  double k[3];
};

// ----------------------------------------------------------------------
// psc_pulse_flattop

struct psc_pulse_flattop {
  double xm, ym, zm; // location of pulse center at time 0 in m 
  double dxm, dym, dzm; // slope of pulse in m
  double zb; // width of pulse in m
  double amplitude_p;   // max amplitude, p-polarization
  double amplitude_s;   // max amplitude, s-polarization
  double phase_p;       // CEP-phase  (from -pi to pi)
  double phase_s;       // CEP-phase  (from -pi to pi)
};

#endif
