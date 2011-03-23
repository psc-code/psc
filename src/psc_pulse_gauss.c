
#include "psc_pulse_private.h"

#include <mrc_params.h>
#include <math.h>

struct psc_pulse_gauss {
  double xm, ym, zm; // location of pulse center at time 0 in m 
  double dxm, dym, dzm; // width of pulse in m
  double amplitude_p;   // max amplitude, p-polarization
  double amplitude_s;   // max amplitude, s-polarization
  double phase_p;       // CEP-phase  (from -pi to pi)
  double phase_s;       // CEP-phase  (from -pi to pi)
  double k[3];
};

// ======================================================================
// psc_pulse_gauss
//

static void
psc_pulse_gauss_setup(struct psc_pulse *pulse)
{
  struct psc_pulse_gauss *gauss = mrc_to_subobj(pulse, struct psc_pulse_gauss);

  // normalization
  gauss->xm /= psc.coeff.ld;
  gauss->ym /= psc.coeff.ld;
  gauss->zm /= psc.coeff.ld;
  gauss->dxm /= psc.coeff.ld;
  gauss->dym /= psc.coeff.ld;
  gauss->dzm /= psc.coeff.ld;
}

static void
psc_pulse_gauss_field(struct psc_pulse *pulse,
		      double xx, double yy, double zz, double tt,
		      double *phase, double *envelope)
{
  struct psc_pulse_gauss *gauss = mrc_to_subobj(pulse, struct psc_pulse_gauss);

  double xr = xx - gauss->xm;
  double yr = yy - gauss->ym;
  double zr = zz - gauss->zm;

  double xl = xr - gauss->k[0] * tt;
  double yl = yr - gauss->k[1] * tt;
  double zl = zr - gauss->k[2] * tt;

  *phase = gauss->k[0] * xr + gauss->k[1] * yr + gauss->k[2] * zr - tt;
  *envelope = (exp(-sqr(xl/gauss->dxm)) *
	       exp(-sqr(yl/gauss->dym)) *
	       exp(-sqr(zl/gauss->dzm)));
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
  { "k"               , VAR(k)            , PARAM_DOUBLE3(0., 0., 0.)       },
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
