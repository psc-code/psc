
#include "psc.h"
#include "psc_case_private.h"
#include "psc_pulse.h"
#include "psc_push_fields.h"
#include "psc_bnd_fields.h"

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

struct psc_case_foils {
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

#define VAR(x) (void *)offsetof(struct psc_case_foils, x)

static struct param psc_case_foils_descr[] = {
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
psc_case_foils_create(struct psc_case *_case)
{
  
  //  float Coeff_FWHM = 0.84932;   // coefficient for putting in values in FWHM of intensity = 1/sqrt(2ln2)
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

#if 0
  struct psc_pulse_gauss prm_s = {
    .xm = 10.   * 1e-6,
    .ym = 2.5   * 1e-6,
    .zm = -45. * 1e-6,
    .dxm = 3.75   * 1e-6 * Coeff_FWHM,
    .dym = 2.   * 1e-6,
    .dzm = 15.   * 1e-6 * Coeff_FWHM,
//    .zb  = 10. * 1e-6,
    .amplitude_s = 1.,
  };
#endif

//  ppsc->pulse_z1 = psc_pulse_flattop_create(&prm_s);
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse_z1 = psc_bnd_fields_get_pulse_z1(bnd_fields);
  psc_pulse_set_type(pulse_z1, "gauss");
}



static void
psc_case_foils_set_from_options(struct psc_case *_case)
{
  ppsc->prm.nmax = 15000;
  ppsc->prm.cpum = 15000;
  ppsc->prm.lw = 1. * 1e-6;
  ppsc->prm.i0 = 2e24;
   
  // n_cr = 1.1e27 for 1 micron wavelength and scales as lambda^-2

  ppsc->prm.n0 = 1.1e29;

  

  ppsc->prm.nicell = 200;

  ppsc->domain.length[0] = 10.0 * 1e-6;			// length of the domain in x-direction (transverse)
  ppsc->domain.length[1] = 0.02 * 1e-6;
  ppsc->domain.length[2] = 10.0  * 1e-6;			// length of the domain in z-direction (longitudinal)

  ppsc->domain.gdims[0] = 100;
  ppsc->domain.gdims[1] = 1;
  ppsc->domain.gdims[2] = 100;

  ppsc->domain.bnd_fld_lo[0] = 1;
  ppsc->domain.bnd_fld_hi[0] = 1;
  ppsc->domain.bnd_fld_lo[1] = 1;
  ppsc->domain.bnd_fld_hi[1] = 1;
  ppsc->domain.bnd_fld_lo[2] = 3; // time
  ppsc->domain.bnd_fld_hi[2] = 2; // upml
  ppsc->domain.bnd_part_lo[0] = 0;
  ppsc->domain.bnd_part_hi[0] = 0;
  ppsc->domain.bnd_part_lo[1] = 0;
  ppsc->domain.bnd_part_hi[1] = 0;
  ppsc->domain.bnd_part_lo[2] = 0;
  ppsc->domain.bnd_part_hi[2] = 0;
}

static void
psc_case_foils_init_field(struct psc_case *_case, mfields_c_t *flds)
{
#if 0
  // FIXME, do we need the ghost points?
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      double dx = ppsc->dx[0], dy = ppsc->dx[1], dz = ppsc->dx[2], dt = ppsc->dt;
      double xx = CRDX(patch, jx), yy = CRDY(patch, jy), zz = CRDZ(patch, jz);
    
      // FIXME, why this time?
      F3(EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, 0.*dt);
      F3(BX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
      F3(EX, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz, 0.*dt);
      F3(BY, jx,jy,jz) = psc_s_pulse_z1(xx + .5*dx, yy, zz + .5*dz, 0.*dt);
    } foreach_3d_g_end;
  }
#endif
}

static void
psc_case_foils_init_npt(struct psc_case *_case, int kind, double x[3], 
			 struct psc_particle_npt *npt)
{
  struct psc_case_foils *foils = mrc_to_subobj(_case, struct psc_case_foils);

  real Te = foils->Te, Ti = foils->Ti;

  real ld = ppsc->coeff.ld;   
 
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

struct psc_case_ops psc_case_foils_ops = {
  .name             = "foils",
  .size             = sizeof(struct psc_case_foils),
  .param_descr      = psc_case_foils_descr,
  .create           = psc_case_foils_create,
  .set_from_options = psc_case_foils_set_from_options,
  .init_field       = psc_case_foils_init_field,
  .init_npt         = psc_case_foils_init_npt,
};
