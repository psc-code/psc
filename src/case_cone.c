
#include "psc.h"
#include "psc_case_private.h"

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

struct psc_case_cone {
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

#define VAR(x) (void *)offsetof(struct psc_case_cone, x)

static struct param psc_case_cone_descr[] = {
  { "Line0_x0"                           , VAR(Line0_x0)                   , PARAM_DOUBLE(11. * 1e-6)           },
  { "Line0_x1"                           , VAR(Line0_x1)                   , PARAM_DOUBLE(7.5 * 1e-6)             },
  { "Line0_z0"                           , VAR(Line0_z0)                   , PARAM_DOUBLE(0.5 * 1e-6)            },
  { "Line0_z1"                           , VAR(Line0_z1)                   , PARAM_DOUBLE(30.5 * 1e-6)           },
  { "Line0_Thickness"                    , VAR(Line0_Thickness)            , PARAM_DOUBLE(0.5 * 1e-6)          },
  { "Line0_Preplasma"                    , VAR(Line0_Preplasma)            , PARAM_DOUBLE(1. * 1e-9)     },
  { "Line1_x0"                           , VAR(Line1_x0)                   , PARAM_DOUBLE(1. * 1e-6)           },
  { "Line1_x1"                           , VAR(Line1_x1)                   , PARAM_DOUBLE(4.5 * 1e-6)             },
  { "Line1_z0"                           , VAR(Line1_z0)                   , PARAM_DOUBLE(0.5 * 1e-6)            },
  { "Line1_z1"                           , VAR(Line1_z1)                   , PARAM_DOUBLE(30.5 * 1e-6)           },
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
psc_case_cone_create(struct psc_case *_case)
{
  
  float Coeff_FWHM = 0.84932;   // coefficient for putting in values in FWHM of intensity = 1/sqrt(2ln2)
   //  float Coeff_FHHM = 1.0;       // uncomment this line if you want the old input  
   // T_L = 3.33 e -15   for 1 micron wavelength
   // T_L = 2.66 e -15   for 0.8 micron wavelength

  struct psc_pulse_gauss prm_p = {
    .xm = 6.   * 1e-6,                       // transverse position of the focus
    .ym = 0.   * 1e-6,
    .zm = -45. * 1e-6,
    .dxm = 3.  * 1e-6 * Coeff_FWHM,
    .dym = 2.   * 1e-6,
    .dzm = 7.5   * 1e-6 * Coeff_FWHM,
//    .zb  = 10. * 1e-6,
    .phase = 0.0,
  };

//  psc.pulse_p_z1 = psc_pulse_flattop_create(&prm_p);
//  psc.pulse_p_z1 = psc_pulse_gauss_create(&prm_p);

#if 0
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
#endif


//  psc.pulse_s_z1 = psc_pulse_flattop_create(&prm_s);
  psc.pulse_p_z1 = psc_pulse_gauss_create(&prm_p);
}



static void
psc_case_cone_set_from_options(struct psc_case *_case)
{
  psc.prm.nmax = 15000;
  psc.prm.cpum = 15000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 2e24;
   
  // n_cr = 1.1e27 for 1 micron wavelength and scales as lambda^-2

  psc.prm.n0 = 1.1e29;

  

  psc.prm.nicell = 200;

  psc.domain.length[0] = 12.0 * 1e-6;			// length of the domain in x-direction (transverse)
  psc.domain.length[1] = 0.02 * 1e-6;
  psc.domain.length[2] = 31.0  * 1e-6;			// length of the domain in z-direction (longitudinal)

  psc.domain.gdims[0] = 1200;
  psc.domain.gdims[1] = 1;
  psc.domain.gdims[2] = 3100;

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
psc_case_cone_init_field(struct psc_case *_case, mfields_base_t *flds)
{
#if 0
  // FIXME, do we need the ghost points?
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2], dt = psc.dt;
      double xx = CRDX(patch, xx), yy = CRDY(patch, yy), zz = CRDZ(patch, zz);
      
      // FIXME, why this time?
      F3_BASE(pf, EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, 0.*dt);
      F3_BASE(pf, BX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
      F3_BASE(pf, EX, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz, 0.*dt);
      F3_BASE(pf, BY, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz + .5*dz, 0.*dt);
    } foreach_3d_g_end;
  }
#endif
}

static void
psc_case_cone_init_npt(struct psc_case *_case, int kind, double x[3], 
			struct psc_particle_npt *npt)
{
  struct psc_case_cone *cone = mrc_to_subobj(_case, struct psc_case_cone);

  real Te = cone->Te, Ti = cone->Ti;

  real ld = psc.coeff.ld;   
 
  real Line0_x0 = cone->Line0_x0 / ld;  
  real Line0_x1 = cone->Line0_x1 / ld;   
  real Line0_z0 = cone->Line0_z0 / ld;
  real Line0_z1 = cone->Line0_z1 / ld;
  real Line0_Thickness = cone->Line0_Thickness / ld;
  real Line0_Preplasma = cone->Line0_Preplasma / ld;

  real Line1_x0 = cone->Line1_x0 / ld;  
  real Line1_x1 = cone->Line1_x1 / ld;   
  real Line1_z0 = cone->Line1_z0 / ld;
  real Line1_z1 = cone->Line1_z1 / ld;
  real Line1_Thickness = cone->Line1_Thickness / ld;
  real Line1_Preplasma = cone->Line1_Preplasma / ld;

  real dens = Line_dens(Line0_x0, Line0_z0, Line0_x1, Line0_z1, x[0], x[2], Line0_Thickness, Line0_Preplasma);
  dens += Line_dens(Line1_x0, Line1_z0, Line1_x1, Line1_z1, x[0], x[2], Line1_Thickness, Line1_Preplasma);

#if 0
  real HollowSphere0_x0 = cone->HollowSphere0_x0 / ld;
  real HollowSphere0_y0 = cone->HollowSphere0_y0 / ld;
  real HollowSphere0_z0 = cone->HollowSphere0_y0 / ld; 
  real HollowSphere0_Radius = cone->R_curv0 / ld;
  real HollowSphere0_Preplasma = cone->HollowSphere0_Preplasma / ld;
  real HollowSphere0_Thickness = cone->HollowSphere0_Thickness / ld;
   
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
    npt->q = cone->charge_state;
    npt->m = cone->mass_ratio;
    npt->n = dens/cone->charge_state;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_cone_ops = {
  .name             = "cone",
  .size             = sizeof(struct psc_case_cone),
  .param_descr      = psc_case_cone_descr,
  .create           = psc_case_cone_create,
  .set_from_options = psc_case_cone_set_from_options,
  .init_field       = psc_case_cone_init_field,
  .init_npt         = psc_case_cone_init_npt,
};
