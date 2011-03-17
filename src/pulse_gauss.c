
#include "psc.h"
#include <mrc_params.h>

#include <math.h>

// ======================================================================
// psc_p_pulse_z1
//
// Laser pulse initialization (p-polarization)
//
// NOTE: The pulse is placed behind of the
// simulation box at a distance "zm" from the
// origin. The pulse then propagates into the 
// simulation box from the left. 
//
//
//  COORDINATE SYSTEM
//
//                          zm        ^ y
//                 <----------------->|
//                                    |
//            laser pulse             |
//                                    |     simulation
//               | | |                |     box
//               | | |   ----->   ^   |
//               | | |         ym |   |
//                                |   |
//          ------------------------------------------------->
//                              (i1n,i2n,i3n)=box origin    z 

#define VAR(x) (void *)offsetof(struct psc_pulse_gauss, x)

static struct param psc_pulse_gauss_descr[] = {
  { "pulse_xm"      , VAR(xm)              , PARAM_DOUBLE(5.  * 1e-6)     },
  { "pulse_ym"      , VAR(ym)              , PARAM_DOUBLE(5.  * 1e-6)     },
  { "pulse_zm"      , VAR(zm)              , PARAM_DOUBLE(-4. * 1e-6)     },
  { "pulse_dxm"     , VAR(dxm)             , PARAM_DOUBLE(1.5 * 1e-6)     },
  { "pulse_dym"     , VAR(dym)             , PARAM_DOUBLE(1.5 * 1e-6)     },
  { "pulse_dzm"     , VAR(dzm)             , PARAM_DOUBLE(1.5 * 1e-6)     },
  { "pulse_phase"   , VAR(phase)           , PARAM_DOUBLE(0.0)            },
  { "pulse_kx"      , VAR(k[0])            , PARAM_DOUBLE(0.)             },
  { "pulse_ky"      , VAR(k[1])            , PARAM_DOUBLE(0.)             },
  { "pulse_kz"      , VAR(k[2])            , PARAM_DOUBLE(0.)             },
  {},
};

#undef VAR

static void
psc_pulse_gauss_setup(struct psc_pulse *pulse)
{
  struct psc_pulse_gauss *prm = pulse->ctx;

  // normalization
  prm->xm /= psc.coeff.ld;
  prm->ym /= psc.coeff.ld;
  prm->zm /= psc.coeff.ld;
  prm->dxm /= psc.coeff.ld;
  prm->dym /= psc.coeff.ld;
  prm->dzm /= psc.coeff.ld;
}

static double
psc_pulse_gauss_field(struct psc_pulse *pulse,
		      double xx, double yy, double zz, double tt)
{
  struct psc_pulse_gauss *prm = pulse->ctx;

  double xr = xx - prm->xm;
  double yr = yy - prm->ym;
  double zr = zz - prm->zm;

  double xl = xr - prm->k[0] * tt;
  double yl = yr - prm->k[1] * tt;
  double zl = zr - prm->k[2] * tt;

  return sin(prm->k[0] * xr + prm->k[1] * yr + prm->k[2] * zr - tt)
    * exp(-sqr(xl/prm->dxm))
    * exp(-sqr(yl/prm->dym))
    * exp(-sqr(zl/prm->dzm));
}

static struct psc_pulse_ops psc_pulse_ops_gauss = {
  .name       = "p_z1_short",
  .ctx_size   = sizeof(struct psc_pulse_gauss),
  .ctx_descr  = psc_pulse_gauss_descr,
  .setup      = psc_pulse_gauss_setup,
  .field      = psc_pulse_gauss_field,
};


struct psc_pulse *
psc_pulse_gauss_create(struct psc_pulse_gauss *ctx)
{
  return psc_pulse_create(&psc_pulse_ops_gauss, ctx);
}

