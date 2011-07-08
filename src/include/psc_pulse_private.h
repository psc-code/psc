
#ifndef PSC_PULSE_PRIVATE_H
#define PSC_PULSE_PRIVATE_H

#include "psc_pulse.h"

struct psc_pulse {
  struct mrc_obj obj;
  double k[3];
  bool is_setup;
};

///Interface for custom pulses
struct psc_pulse_ops {
  MRC_SUBCLASS_OPS(struct psc_pulse);
  double (*field_s)(struct psc_pulse *,
		    double x, double y, double z, double t);	///< The field value in s-polarisation
  double (*field_p)(struct psc_pulse *,
		    double x, double y, double z, double t);	///< The field value in p-polarisation
};  

extern struct psc_pulse_ops psc_pulse_none_ops;
extern struct psc_pulse_ops psc_pulse_gauss_ops;
extern struct psc_pulse_ops psc_pulse_flattop_ops;

#define psc_pulse_ops(pulse) ((struct psc_pulse_ops *)((pulse)->obj.ops))

#endif

