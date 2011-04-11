
#include "psc_pulse_private.h"

#include <math.h>

struct psc_pulse_flattop {
  double xm[3];         // location of pulse center at time 0 in m 
  double dxm[3];        // slope of pulse in m
  double zb;            // width of pulse in m
  double amplitude_p;   // max amplitude, p-polarization
  double amplitude_s;   // max amplitude, s-polarization
  double phase_p;       // CEP-phase  (from -pi to pi)
  double phase_s;       // CEP-phase  (from -pi to pi)
};

static void
psc_pulse_flattop_setup(struct psc_pulse *pulse)
{
  struct psc_pulse_flattop *flattop = mrc_to_subobj(pulse, struct psc_pulse_flattop);

  // normalization
  for (int d = 0; d < 3; d++) {
    flattop->xm[d] /= psc.coeff.ld;
    flattop->dxm[d] /= psc.coeff.ld;
  }
  flattop->zb /= psc.coeff.ld;
}

void
psc_pulse_flattop_field(struct psc_pulse *pulse, 
			double xx, double yy, double zz, double tt,
			double *phase, double *envelope)
{
  struct psc_pulse_flattop *flattop = mrc_to_subobj(pulse, struct psc_pulse_flattop);

  // FIXME, only for z direction right now

  double xl = xx;
  double yl = yy;
  double zl = zz - tt;

  double xr = xl - flattop->xm[0];
  double yr = yl - flattop->xm[1];
  double zr = zl - flattop->xm[2];

  *phase = zr;
  *envelope = (exp(-sqr(xr/flattop->dxm[0])) *
	       exp(-sqr(yr/flattop->dxm[1])) *
	       1. / (1.+exp((fabs(zr)-flattop->zb)/flattop->dxm[2])));
}

static double
psc_pulse_flattop_field_p(struct psc_pulse *pulse,
			  double xx, double yy, double zz, double tt)
{
  struct psc_pulse_flattop *flattop = mrc_to_subobj(pulse, struct psc_pulse_flattop);

  double phase, envelope;
  psc_pulse_flattop_field(pulse, xx, yy, zz, tt, &phase, &envelope);

  return flattop->amplitude_p * envelope * sin(phase + flattop->phase_p);
}

static double
psc_pulse_flattop_field_s(struct psc_pulse *pulse,
			  double xx, double yy, double zz, double tt)
{
  struct psc_pulse_flattop *flattop = mrc_to_subobj(pulse, struct psc_pulse_flattop);

  double phase, envelope;
  psc_pulse_flattop_field(pulse, xx, yy, zz, tt, &phase, &envelope);

  return flattop->amplitude_s * envelope * sin(phase + flattop->phase_s);
}

#define VAR(x) (void *)offsetof(struct psc_pulse_flattop, x)
static struct param psc_pulse_flattop_descr[] = {
  { "m"               , VAR(xm)           , PARAM_DOUBLE3(0., 0., 0.)       },
  { "dm"              , VAR(dxm)          , PARAM_DOUBLE3(1e-6, 1e-6, 1e-6) },
  { "zb"              , VAR(zb)           , PARAM_DOUBLE(1.5 * 1e-6)        },
  { "phase_p"         , VAR(phase_p)      , PARAM_DOUBLE(0.0)               },
  { "phase_s"         , VAR(phase_s)      , PARAM_DOUBLE(0.0)               },
  { "amplitude_p"     , VAR(amplitude_p)  , PARAM_DOUBLE(0.)                },
  { "amplitude_s"     , VAR(amplitude_s)  , PARAM_DOUBLE(0.)                },
  {},
};
#undef VAR

struct psc_pulse_ops psc_pulse_flattop_ops = {
  .name        = "flattop",
  .size        = sizeof(struct psc_pulse_flattop),
  .param_descr = psc_pulse_flattop_descr,
  .setup       = psc_pulse_flattop_setup,
  .field_p     = psc_pulse_flattop_field_p,
  .field_s     = psc_pulse_flattop_field_s,
};
