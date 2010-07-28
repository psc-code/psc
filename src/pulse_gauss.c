
#include "psc.h"
#include "util/params.h"

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

#define VAR(x) (void *)offsetof(struct psc_p_pulse_z1_param, x)

static struct param psc_p_pulse_z1_descr[] = {
  { "pulse_xm"      , VAR(xm)              , PARAM_DOUBLE(5.  * 1e-6)     },
  { "pulse_ym"      , VAR(ym)              , PARAM_DOUBLE(5.  * 1e-6)     },
  { "pulse_zm"      , VAR(zm)              , PARAM_DOUBLE(-4. * 1e-6)     },
  { "pulse_dxm"     , VAR(dxm)             , PARAM_DOUBLE(1.5 * 1e-6)     },
  { "pulse_dym"     , VAR(dym)             , PARAM_DOUBLE(1.5 * 1e-6)     },
  { "pulse_dzm"     , VAR(dzm)             , PARAM_DOUBLE(1.5 * 1e-6)     },
  {},
};

#undef VAR

static void
psc_pulse_p_z1_setup(struct psc_pulse *pulse)
{
  struct psc_p_pulse_z1_param *prm = pulse->ctx;

  // normalization
  prm->xm /= psc.coeff.ld;
  prm->ym /= psc.coeff.ld;
  prm->zm /= psc.coeff.ld;
  prm->dxm /= psc.coeff.ld;
  prm->dym /= psc.coeff.ld;
  prm->dzm /= psc.coeff.ld;

  pulse->is_setup = true;
}

static double
psc_pulse_p_z1_short_field(struct psc_pulse *pulse,
			   double xx, double yy, double zz, double tt)
{
  struct psc_p_pulse_z1_param *prm = pulse->ctx;

  if (!pulse->is_setup) {
    psc_pulse_p_z1_setup(pulse);
  }

  //  double xl = xx;
  double yl = yy;
  double zl = zz - tt;

  //  double xr = xl - prm->xm;
  double yr = yl - prm->ym;
  double zr = zl - prm->zm;

  return sin(zr)
    // * exp(-sqr(xr/prm->dxm))
    * exp(-sqr(yr/prm->dym))
    * exp(-sqr(zr/prm->dzm));
}

static struct psc_pulse_ops psc_pulse_ops_p_z1_short = {
  .name       = "p_z1_short",
  .ctx_size   = sizeof(struct psc_p_pulse_z1_param),
  .field      = psc_pulse_p_z1_short_field,
};


struct psc_pulse *
psc_pulse_p_z1_short_create(struct psc_p_pulse_z1_param *prm)
{
  struct psc_pulse *pulse = psc_pulse_create(&psc_pulse_ops_p_z1_short);

  struct psc_p_pulse_z1_param *ctx = pulse->ctx;
  if (prm) { // custom defaults were passed
    memcpy(ctx, prm, sizeof(*ctx));
    params_parse_cmdline_nodefault(ctx, psc_p_pulse_z1_descr,
				   "PSC P pulse z1", MPI_COMM_WORLD);
  } else {
    params_parse_cmdline(ctx, psc_p_pulse_z1_descr,
			 "PSC P pulse z1", MPI_COMM_WORLD);
  }
  params_print(ctx, psc_p_pulse_z1_descr, "PSC P pulse z1",
	       MPI_COMM_WORLD);

  return pulse;
}

