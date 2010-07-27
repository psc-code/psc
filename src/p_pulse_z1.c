
#include "psc.h"

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

double
p_pulse_z1(double xx, double yy, double zz, double tt)
{
  double dxm = 5.0*1.0e-6;
  double dym = 5.0*1.0e-6;
  double dzm = 1.0*1.0e-6;
  double xm = 2.0*1.0e-5;
  double ym = 2.0*1.0e-5;
  double zm =-2.0*1.0e-6;

  // normalization
  xm /= psc.coeff.ld;
  ym /= psc.coeff.ld;
  zm /= psc.coeff.ld;
  dxm /= psc.coeff.ld;
  dym /= psc.coeff.ld;
  dzm /= psc.coeff.ld;

  //  double xl = xx;
  double yl = yy;
  double zl = zz - tt;

  //  double xr = xl - xm;
  double yr = yl - ym;
  double zr = zl - zm;

  return sin(zr)
    // * exp(-sqr(xr/dxm))
    * exp(-sqr(yr/dym))
    * exp(-sqr(zr/dzm));
}

struct psc_pulse_ops psc_pulse_ops_p_z1_short = {
  .name       = "p_z1_short",
  .p_pulse_z1 = p_pulse_z1,
};

double
p_pulse_flattop_z1(double xx, double yy, double zz, double tt)
{
  double dxm = 4.0*1e-6;
  double dym = 4.0*1e-6;
  double dzm = .1 *1e-6;
  double xm = .01 *1e-6;
  double ym = .01 *1e-6;
  double zm =-301.*1e-6;
  double zb = 300.*1e-6;

  // normalization
  xm /= psc.coeff.ld;
  ym /= psc.coeff.ld;
  zm /= psc.coeff.ld;
  dxm /= psc.coeff.ld;
  dym /= psc.coeff.ld;
  dzm /= psc.coeff.ld;
  zb  /= psc.coeff.ld;

  double xl = xx;
  double yl = yy;
  double zl = zz - tt;

  double xr = xl - xm;
  double yr = yl - ym;
  double zr = zl - zm;

  return sin(zr)
    * exp(-sqr(xr/dxm))
    * exp(-sqr(yr/dym))
    * 1. / (1.+exp((fabs(zr)-zb)/dzm));
}

struct psc_pulse_ops psc_pulse_ops_p_z1_flattop = {
  .name       = "p_z1_flattop",
  .p_pulse_z1 = p_pulse_flattop_z1,
};
