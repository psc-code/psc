
#include "psc.h"
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================
// Foils : different density distribution examples:
// 
// 1. Line
// 2. Circle
// 3. Part of circle
//-----------------------------------------------------
//  
// 1. Line
// needs the coordinates of its beginning and end, thickness 
//

struct foils {
  
  double Line0_x0, Line0_z0;   // coordinates of the beginning of the line0
  double Line0_x1, Line0_z1;   // coordinates of the end of the line0
  double Line0_Thickness;       // thickness of the line
  double Line0_Preplasma;    
 
  double Line1_x0, Line1_z0;
  double Line1_x1, Line1_z1;
  double Line1_Thickness;
  double Line1_Preplasma;

  double Te, Ti;
  
  double HollowSphere0_x0, HollowSphere0_y0, HollowSphere0_z0; // location of density center in m of the first(0) foil
  double HollowSphere0_Preplasma; // gradient of density profile in m of the first foil
  double HollowSphere0_Thickness; // width of transverse / longitudinal 
                                 // density profile in m of the first foil 
  double mass_ratio; // M_i / M_e
  double charge_state;   // Charge state of the ion
  double R_curv0;  // curvature of the hollowsphere0 foil  in meters
};

#define VAR(x) (void *)offsetof(struct foils, x)

static struct param foils_descr[] = {
  { "Line0_x0"                           , VAR(Line0_x0)                   , PARAM_DOUBLE(15. * 1e-6)           },
  { "Line0_x1"                           , VAR(Line0_x1)                   , PARAM_DOUBLE(11.5 * 1e-6)             },
  { "Line0_z0"                           , VAR(Line0_z0)                   , PARAM_DOUBLE(3.0 * 1e-6)            },
  { "Line0_z1"                           , VAR(Line0_z1)                   , PARAM_DOUBLE(23. * 1e-6)           },
  { "Line0_Thickness"                    , VAR(Line0_Thickness)            , PARAM_DOUBLE(0.5 * 1e-6)          },
  { "Line0_Preplasma"                    , VAR(Line0_Preplasma)            , PARAM_DOUBLE(1. * 1e-9)     },
  { "Line1_x0"                           , VAR(Line1_x0)                   , PARAM_DOUBLE(5. * 1e-6)           },
  { "Line1_x1"                           , VAR(Line1_x1)                   , PARAM_DOUBLE(8.5 * 1e-6)             },
  { "Line1_z0"                           , VAR(Line1_z0)                   , PARAM_DOUBLE(3.0 * 1e-6)            },
  { "Line1_z1"                           , VAR(Line1_z1)                   , PARAM_DOUBLE(23. * 1e-6)           },
  { "Line1_Thickness"                    , VAR(Line1_Thickness)            , PARAM_DOUBLE(0.5 * 1e-6)          },
  { "Line1_Preplasma"                    , VAR(Line1_Preplasma)            , PARAM_DOUBLE(1. * 1e-9)     },
  { "Te"                                 , VAR(Te)                         , PARAM_DOUBLE(0.)             },
  { "Ti"                                 , VAR(Ti)                         , PARAM_DOUBLE(0.)             },
  { "HollowSphere0_x0"                   , VAR(HollowSphere0_x0)           , PARAM_DOUBLE(2.5 * 1e-6)     },
  { "HollowSphere0_y0"                   , VAR(HollowSphere0_y0)           , PARAM_DOUBLE(.01 * 1e-6)     },
  { "HollowSphere0_z0"                   , VAR(HollowSphere0_z0)           , PARAM_DOUBLE(2.5  * 1e-6)     },
  { "HollowSphere0_Preplasma"            , VAR(HollowSphere0_Preplasma)    , PARAM_DOUBLE(10.  * 1e-9)     },
  { "HollowSphere0_Thickness"            , VAR(HollowSphere0_Thickness)    , PARAM_DOUBLE(200. * 1e-9)     },
  { "mass_ratio"                         , VAR(mass_ratio)                 , PARAM_DOUBLE(12.*1836.)      },
  { "charge_state"                       , VAR(charge_state)               , PARAM_DOUBLE(6.)             },
  { "R_curv0"                            , VAR(R_curv0)                    , PARAM_DOUBLE(2.5 * 1e-6)             },
  {},
};

#undef VAR

static real Line_dens(double x0, double z0, double x1, double z1, double xc, double zc, double Thickness, double Preplasma)
{
    // returns the density in the current cell for the line density distribution
    // x0,z0 - coordinates of the beginning of the line
    // x1,z1 - coordinates of the end of the line
    // xc, zc - current coordinates of the grid
    // Thickness - thickness of the line
    // Preplasma - preplasma of the line

    real Length = sqrt(sqr(x0-x1)+sqr(z0-z1));
    real cosrot = 1.0;
    real sinrot = 0.0;

    if(Length!=0.0)
    { 
      cosrot = -(z0 - z1) / Length;
      sinrot = -(x0 - x1) / Length;  
    }

  real xmiddle = (x0+x1)*0.5;
  real zmiddle = (z0+z1)*0.5;

  real xr = sinrot * (xc - xmiddle) + cosrot * (zc - zmiddle) + xmiddle;
  real zr = sinrot * (zc - zmiddle) - cosrot * (xc - xmiddle) + zmiddle;

   real argx = (fabs(xr-xmiddle)-0.5*Length)/Preplasma;
   real argz = (fabs(zr-zmiddle)-0.5*Thickness)/Preplasma;
  if (argx > 200.0) argx = 200.0;
  if (argz > 200.0) argz = 200.0;

  return 1. / ((1. + exp(argx)) * (1. + exp(argz)));
}

#if 0
static real HollowSphere_dens(double x0, double z0, double Radius, double xc, double zc, double Thickness, double Preplasma)
{
    // returns the density in the current cell for the hollow sphere density distribution
    // x0,z0 - coordinates of the center of the sphere
    // Radius - radius of the sphere
    // xc, zc - current coordinates of the grid
    // Thickness - thickness of the wall of the sphere
    // Preplasma - preplasma 
   
  real RadiusAtCurrentCell = sqrt((xc-x0)*(xc-x0)+(zc-z0)*(zc-z0));
  real argsphere = (fabs(RadiusAtCurrentCell-Radius)-Thickness)/Preplasma;

  if (argsphere > 200.0) argsphere = 200.0;

  return 1./(1.+exp(argsphere));
}
#endif

static void
foils_create(struct psc_case *Case)
{
  
  float Coeff_FWHM = 0.84932;   // coefficient for putting in values in FWHM of intensity = 1/sqrt(2ln2)
   //  float Coeff_FHHM = 1.0;       // uncomment this line if you want the old input  
   // T_L = 3.33 e -15   for 1 micron wavelength
   // T_L = 2.66 e -15   for 0.8 micron wavelength

#if 0
  struct psc_pulse_gauss prm_p = {
    .xm = 10.   * 1e-6,                       // transverse position of the focus
    .ym = 2.5   * 1e-6,
    .zm = -45. * 1e-6,
    .dxm = 3.75  * 1e-6 * Coeff_FWHM,
    .dym = 2.   * 1e-6,
    .dzm = 15.   * 1e-6 * Coeff_FWHM,
//    .zb  = 10. * 1e-6,
    .phase = 0.0,
  };

//  psc.pulse_p_z1 = psc_pulse_flattop_create(&prm_p);
//  psc.pulse_p_z1 = psc_pulse_gauss_create(&prm_p);
#endif

  struct psc_pulse_gauss prm_s = {
    .xm = 10.   * 1e-6,
    .ym = 2.5   * 1e-6,
    .zm = -45. * 1e-6,
    .dxm = 3.75   * 1e-6 * Coeff_FWHM,
    .dym = 2.   * 1e-6,
    .dzm = 15.   * 1e-6 * Coeff_FWHM,
//    .zb  = 10. * 1e-6,
    .phase = 0.0,
  };

//  psc.pulse_s_z1 = psc_pulse_flattop_create(&prm_s);
  psc.pulse_s_z1 = psc_pulse_gauss_create(&prm_s);
}



static void
foils_init_param(struct psc_case *Case)
{
  psc.prm.nmax = 15000;
  psc.prm.cpum = 15000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 2e24;
   
  // n_cr = 1.1e27 for 1 micron wavelength and scales as lambda^-2

  psc.prm.n0 = 1.1e29;

  

  psc.prm.nicell = 200;

  psc.domain.length[0] = 10.0 * 1e-6;			// length of the domain in x-direction (transverse)
  psc.domain.length[1] = 0.02 * 1e-6;
  psc.domain.length[2] = 10.0  * 1e-6;			// length of the domain in z-direction (longitudinal)

  psc.domain.itot[0] = 100;				// total number of steps in x-direction. dx=length/itot;
  psc.domain.itot[1] = 10;				
  psc.domain.itot[2] = 100;				// total number of steps in z-direction. dz=length/itot;
  psc.domain.ilo[0] = 0;
  psc.domain.ilo[1] = 9;
  psc.domain.ilo[2] = 0;
  psc.domain.ihi[0] = 100;
  psc.domain.ihi[1] = 10;
  psc.domain.ihi[2] = 100;

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
foils_init_field(struct psc_case *Case)
{
#if 0
  // FIXME, do we need the ghost points?
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2], dt = psc.dt;
	double xx = jx * dx, yy = jy * dy, zz = jz * dz;

	// FIXME, why this time?
	FF3(EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, 0.*dt);
	FF3(BX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
	FF3(EX, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz, 0.*dt);
	FF3(BY, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz + .5*dz, 0.*dt);
      }
    }
  }
#endif
}

static void
foils_init_npt(struct psc_case *Case, int kind, double x[3], 
		  struct psc_particle_npt *npt)
{
  struct foils *foils = Case->ctx;

  real Te = foils->Te, Ti = foils->Ti;

  real ld = psc.coeff.ld;   
 
  real Line0_x0 = foils->Line0_x0 / ld;  
  real Line0_x1 = foils->Line0_x1 / ld;   
  real Line0_z0 = foils->Line0_z0 / ld;
  real Line0_z1 = foils->Line0_z1 / ld;
  real Line0_Thickness = foils->Line0_Thickness / ld;
  real Line0_Preplasma = foils->Line0_Preplasma / ld;

  real Line1_x0 = foils->Line1_x0 / ld;  
  real Line1_x1 = foils->Line1_x1 / ld;   
  real Line1_z0 = foils->Line1_z0 / ld;
  real Line1_z1 = foils->Line1_z1 / ld;
  real Line1_Thickness = foils->Line1_Thickness / ld;
  real Line1_Preplasma = foils->Line1_Preplasma / ld;

  real dens = Line_dens(Line0_x0, Line0_z0, Line0_x1, Line0_z1, x[0], x[2], Line0_Thickness, Line0_Preplasma);
  dens += Line_dens(Line1_x0, Line1_z0, Line1_x1, Line1_z1, x[0], x[2], Line1_Thickness, Line1_Preplasma);

#if 0
  real HollowSphere0_x0 = foils->HollowSphere0_x0 / ld;
  real HollowSphere0_y0 = foils->HollowSphere0_y0 / ld;
  real HollowSphere0_z0 = foils->HollowSphere0_y0 / ld; 
  real HollowSphere0_Radius = foils->R_curv0 / ld;
  real HollowSphere0_Preplasma = foils->HollowSphere0_Preplasma / ld;
  real HollowSphere0_Thickness = foils->HollowSphere0_Thickness / ld;
   
  dens += HollowSphere_dens(HollowSphere0_x0, HollowSphere0_z0, HollowSphere0_Radius, x[0], x[2], HollowSphere0_Thickness, HollowSphere0_Preplasma);
#endif

  if (dens>1.0) dens=1.0;

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
    npt->q = foils->charge_state;
    npt->m = foils->mass_ratio;
    npt->n = dens/foils->charge_state;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_ops_foils = {
  .name       = "foils",
  .ctx_size   = sizeof(struct foils),
  .ctx_descr  = foils_descr,
  .create     = foils_create,
  .init_param = foils_init_param,
  .init_field = foils_init_field,
  .init_npt   = foils_init_npt,
};
