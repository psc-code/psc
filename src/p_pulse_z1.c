
#include "psc.h"
#include "util/params.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

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

struct psc_p_pulse_z1 {
  struct psc_pulse pulse;

  double xm, ym, zm; // location of pulse center at time 0 in m 
  double dxm, dym, dzm; // width of pulse in m

  bool is_setup;
};

#define VAR(x) (void *)offsetof(struct psc_p_pulse_z1, x)

static struct param psc_p_pulse_z1_descr[] = {
  { "pulse_xm"      , VAR(xm)              , PARAM_DOUBLE(20. * 1e-6)     },
  { "pulse_ym"      , VAR(ym)              , PARAM_DOUBLE(20. * 1e-6)     },
  { "pulse_zm"      , VAR(zm)              , PARAM_DOUBLE(-2. * 1e-6)     },
  { "pulse_dxm"     , VAR(dxm)             , PARAM_DOUBLE(5.  * 1e-6)     },
  { "pulse_dym"     , VAR(dym)             , PARAM_DOUBLE(5.  * 1e-6)     },
  { "pulse_dzm"     , VAR(dzm)             , PARAM_DOUBLE(1.  * 1e-6)     },
  {},
};

#undef VAR

static void
psc_pulse_p_z1_short_destroy(struct psc_pulse *__pulse)
{
  struct psc_p_pulse_z1 *pulse = (struct psc_p_pulse_z1 *) __pulse;
  free(pulse);
}

static void
psc_pulse_p_z1_setup(struct psc_p_pulse_z1 *pulse)
{
  // normalization
  pulse->xm /= psc.coeff.ld;
  pulse->ym /= psc.coeff.ld;
  pulse->zm /= psc.coeff.ld;
  pulse->dxm /= psc.coeff.ld;
  pulse->dym /= psc.coeff.ld;
  pulse->dzm /= psc.coeff.ld;

  pulse->is_setup = true;
}

static double
psc_pulse_p_z1_short_field(struct psc_pulse *__pulse,
			   double xx, double yy, double zz, double tt)
{
  struct psc_p_pulse_z1 *pulse = (struct psc_p_pulse_z1 *) __pulse;

  if (!pulse->is_setup) {
    psc_pulse_p_z1_setup(pulse);
  }

  //  double xl = xx;
  double yl = yy;
  double zl = zz - tt;

  //  double xr = xl - pulse->xm;
  double yr = yl - pulse->ym;
  double zr = zl - pulse->zm;

  return sin(zr)
    // * exp(-sqr(xr/pulse->dxm))
    * exp(-sqr(yr/pulse->dym))
    * exp(-sqr(zr/pulse->dzm));
}

static struct psc_pulse_ops psc_pulse_ops_p_z1_short = {
  .name       = "p_z1_short",
  .destroy    = psc_pulse_p_z1_short_destroy,
  .field      = psc_pulse_p_z1_short_field,
};

struct psc_pulse *
psc_pulse_p_z1_short_create(void)
{
  struct psc_p_pulse_z1 *pulse = malloc(sizeof(*pulse));
  memset(pulse, 0, sizeof(*pulse));
  pulse->pulse.ops = &psc_pulse_ops_p_z1_short;

  params_parse_cmdline(pulse, psc_p_pulse_z1_descr, "PSC P pulse z1", MPI_COMM_WORLD);
  params_print(pulse, psc_p_pulse_z1_descr, "PSC P pulse z1", MPI_COMM_WORLD);

  return (struct psc_pulse *) pulse;
}

