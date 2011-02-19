
#include "psc.h"
#include "util/params.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================
// curved thinfoil
//
// Thin overdense foil(s) + curvature + laser pulse
//
//              **        in this case the R_curv0 is negative, if R_curv0
//               **       is positive then the foil will be curved to the
//                **      other side
//       R_curv0   **
//        *--------**<-- x0,y0,z0 
//                 **
//                **
//               **
//              **

struct curvedfoil {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m of the first(0) foil
  double L0; // gradient of density profile in m of the first foil
  double width0; // width of transverse / longitudinal 
                                 // density profile in m of the first foil 
  double mass_ratio; // M_i / M_e
  double R_curv0;  // curvature of the first foil  in meters
};

#define VAR(x) (void *)offsetof(struct curvedfoil, x)

static struct param curvedfoil_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.)             },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.)             },
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(2.5 * 1e-6)     },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(.01 * 1e-6)     },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(2.5  * 1e-6)     },
  { "L0"            , VAR(L0)              , PARAM_DOUBLE(10.  * 1e-9)     },
  { "width0"        , VAR(width0)          , PARAM_DOUBLE(200. * 1e-9)     },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(12.*1836.)      },
  { "R_curv0"     , VAR(R_curv0)      ,       PARAM_DOUBLE(2.5 * 1e-6)             },
  {},
};

#undef VAR

static void
curvedfoil_create(struct psc_case *Case)
{
  struct psc_pulse_gauss prm_p = {
    .xm = 2.5   * 1e-6,
    .ym = 2.5   * 1e-6,
    .zm = -0. * 1e-6,
    .dxm = 0.5   * 1e-6,
    .dym = 2.   * 1e-6,
    .dzm = 4.   * 1e-6,
//    .zb  = 10. * 1e-6,
    .phase = 0.0,
  };
  
  struct psc_pulse_gauss prm_s = {
    .xm = 2.5   * 1e-6,
    .ym = 2.5   * 1e-6,
    .zm = -0. * 1e-6,
    .dxm = 0.5   * 1e-6,
    .dym = 2.   * 1e-6,
    .dzm = 4.   * 1e-6,
//    .zb  = 10. * 1e-6,
    .phase = M_PI / 2.,
  };

//  psc.pulse_p_z1 = psc_pulse_flattop_create(&prm_p);
//  psc.pulse_s_z1 = psc_pulse_flattop_create(&prm_s);

  psc.pulse_p_z1 = psc_pulse_gauss_create(&prm_p);
  psc.pulse_s_z1 = psc_pulse_gauss_create(&prm_s);
}



static void
curvedfoil_init_param(struct psc_case *Case)
{
  psc.prm.nmax = 500;
  psc.prm.cpum = 25000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 1.0e20;
  psc.prm.n0 = 2.0e29;

  psc.prm.nicell = 200;

  psc.domain.length[0] = 5 * 1e-6;			// length of the domain in x-direction (transverse)
  psc.domain.length[1] = 0.02 * 1e-6;
  psc.domain.length[2] = 5.0  * 1e-6;			// length of the domain in z-direction (longitudinal)

  psc.domain.itot[0] = 200;				// total number of steps in x-direction. dx=length/itot;
  psc.domain.itot[1] = 20;				
  psc.domain.itot[2] = 200;				// total number of steps in z-direction. dz=length/itot;
  psc.domain.ilo[0] = 0;
  psc.domain.ilo[1] = 9;
  psc.domain.ilo[2] = 0;
  psc.domain.ihi[0] = 200;
  psc.domain.ihi[1] = 10;
  psc.domain.ihi[2] = 200;

  psc.domain.bnd_fld_lo[0] = 1;
  psc.domain.bnd_fld_hi[0] = 1;
  psc.domain.bnd_fld_lo[1] = 1;
  psc.domain.bnd_fld_hi[1] = 1;
  psc.domain.bnd_fld_lo[2] = 3; // time
  psc.domain.bnd_fld_hi[2] = 2; // upml
  psc.domain.bnd_part[0] = 0;
  psc.domain.bnd_part[1] = 0;
  psc.domain.bnd_part[2] = 0;
}

static void
curvedfoil_init_field(struct psc_case *Case)
{
  // FIXME, do we need the ghost points?
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2], dt = psc.dt;
	double xx = jx * dx, yy = jy * dy, zz = jz * dz;

	// FIXME, why this time?
	F3_BASE(EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, 0.*dt);
	F3_BASE(BX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
	F3_BASE(EX, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz, 0.*dt);
	F3_BASE(BY, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz + .5*dz, 0.*dt);
      }
    }
  }
}

#if 0
static void
curvedfoil_init_npt(struct psc_case *Case, int kind, double x[3], 
		  struct psc_particle_npt *npt)
{
  struct curvedfoil *curvedfoil = Case->ctx;

  real Te = curvedfoil->Te, Ti = curvedfoil->Ti;

  real ld = psc.coeff.ld;    // what is ld?

  real x0 = curvedfoil->x0 / ld;
  real y0 = curvedfoil->y0 / ld;
  real z0 = curvedfoil->z0 / ld;
  real L0 = curvedfoil->L0 / ld;
  real width0 = curvedfoil->width0 / ld;
  real R_curv0 = curvedfoil->R_curv0 / ld;

  //real yr = x[1];
  //real xr = cos(rot) * (x[0]-x0) - sin(rot) * (x[2]-z0) + x0;
  //real zr = cos(rot) * (x[2]-z0) + sin(rot) * (x[0]-x0) + z0;

  real xr = x[0];
  real yr = x[1];
  real zr = x[2];

  real xsphere0 = x0;                    // coordinates 
  real ysphere0 = y0;                    // of the 
  real zsphere0 = z0;                    // sphere0 center

  real Radius0AtCurrentCell = sqrt((xr-xsphere0)*(xr-xsphere0)+(yr-ysphere0)*(yr-ysphere0)+(zr-zsphere0)*(zr-zsphere0));
  real argsphere0 = (fabs(Radius0AtCurrentCell-R_curv0)-width0)/L0;

  if (argsphere0 > 200.0) argsphere0 = 200.0;

  real dens = 1./(1.+exp(argsphere0));

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->n = dens;
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = curvedfoil->mass_ratio;
    npt->n = dens;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}
#endif

struct psc_case_ops psc_case_ops_curvedfoil = {
  .name       = "curvedfoil",
  .ctx_size   = sizeof(struct curvedfoil),
  .ctx_descr  = curvedfoil_descr,
  .create     = curvedfoil_create,
  .init_param = curvedfoil_init_param,
  .init_field = curvedfoil_init_field,
//  .init_npt   = curvedfoil_init_npt,
};
