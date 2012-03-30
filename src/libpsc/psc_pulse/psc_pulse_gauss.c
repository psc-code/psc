
#include "psc_pulse_private.h"

#include <mrc_params.h>
#include <math.h>

struct psc_pulse_gauss {
  double xm[3];         // location of pulse center at time 0 in m 
  double dxm[3];        // width of pulse in m
  double amplitude_p;   // max amplitude, p-polarization
  double amplitude_s;   // max amplitude, s-polarization
  double phase_p;       // CEP-phase  (from -pi to pi)
  double phase_s;       // CEP-phase  (from -pi to pi)
};

// ======================================================================
// psc_pulse_gauss
//

static void
psc_pulse_gauss_setup(struct psc_pulse *pulse)
{
  psc_pulse_setup_super(pulse);

  // FIXME, not good to change parameters, messes up ::view()
  struct psc_pulse_gauss *gauss = mrc_to_subobj(pulse, struct psc_pulse_gauss);

  // normalization
  for (int d = 0; d < 3; d++) {
    gauss->xm[d] /= ppsc->coeff.ld;
    gauss->dxm[d] /= ppsc->coeff.ld;
  }
}

static void
psc_pulse_gauss_field(struct psc_pulse *pulse,
		      double xx, double yy, double zz, double tt,
		      double *phase, double *envelope)
{
  struct psc_pulse_gauss *gauss = mrc_to_subobj(pulse, struct psc_pulse_gauss);

  double xr = xx - gauss->xm[0];
  double yr = yy - gauss->xm[1];
  double zr = zz - gauss->xm[2];

  double xl = xr - pulse->k[0] * tt;
  double yl = yr - pulse->k[1] * tt;
  double zl = zr - pulse->k[2] * tt;

  *phase = pulse->k[0] * xr + pulse->k[1] * yr + pulse->k[2] * zr - tt;
  *envelope = (exp(-sqr(xl/gauss->dxm[0])) *
	       exp(-sqr(yl/gauss->dxm[1])) *
	       exp(-sqr(zl/gauss->dxm[2])));
}

static double
psc_pulse_gauss_field_p(struct psc_pulse *pulse,
			double xx, double yy, double zz, double tt)
{
  struct psc_pulse_gauss *gauss = mrc_to_subobj(pulse, struct psc_pulse_gauss);

  double phase, envelope;
  psc_pulse_gauss_field(pulse, xx, yy, zz, tt, &phase, &envelope);

  return gauss->amplitude_p * envelope * sin(phase + gauss->phase_p);
}

static double
psc_pulse_gauss_field_s(struct psc_pulse *pulse,
			double xx, double yy, double zz, double tt)
{
  struct psc_pulse_gauss *gauss = mrc_to_subobj(pulse, struct psc_pulse_gauss);

  double phase, envelope;
  psc_pulse_gauss_field(pulse, xx, yy, zz, tt, &phase, &envelope);

  return gauss->amplitude_s * envelope * sin(phase + gauss->phase_s);
}

#define VAR(x) (void *)offsetof(struct psc_pulse_gauss, x)
static struct param psc_pulse_gauss_descr[] = {
  { "m"               , VAR(xm)           , PARAM_DOUBLE3(0., 0., 0.)       },
  { "dm"              , VAR(dxm)          , PARAM_DOUBLE3(1e-6, 1e-6, 1e-6) },
  { "phase_p"         , VAR(phase_p)      , PARAM_DOUBLE(0.0)               },
  { "phase_s"         , VAR(phase_s)      , PARAM_DOUBLE(0.0)               },
  { "amplitude_p"     , VAR(amplitude_p)  , PARAM_DOUBLE(0.)                },
  { "amplitude_s"     , VAR(amplitude_s)  , PARAM_DOUBLE(0.)                },
  {},
};
#undef VAR

struct psc_pulse_ops psc_pulse_gauss_ops = {
  .name        = "gauss",
  .size        = sizeof(struct psc_pulse_gauss),
  .param_descr = psc_pulse_gauss_descr,
  .setup       = psc_pulse_gauss_setup,
  .field_p     = psc_pulse_gauss_field_p,
  .field_s     = psc_pulse_gauss_field_s,
};
