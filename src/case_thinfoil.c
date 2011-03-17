
#include "psc.h"
#include "psc_case_private.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================
// thinfoil
//
// FIXME description

struct psc_case_thinfoil {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double Lx, Ly, Lz; // gradient of density profile in m
  double widthx, widthy, widthz; // width of transverse / longitudinal 
                                 // density profile in m
  double rot;        // rotation angle in degrees
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct psc_case_thinfoil, x)

static struct param psc_case_thinfoil_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.1)            },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.1)            },
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(.01 * 1e-6)     },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(.01 * 1e-6)     },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(1.  * 1e-6)     },
  { "Lx"            , VAR(Lx)              , PARAM_DOUBLE(1.  * 1e-9)     },
  { "Ly"            , VAR(Ly)              , PARAM_DOUBLE(1.  * 1e-9)     },
  { "Lz"            , VAR(Lz)              , PARAM_DOUBLE(.5  * 1e-9)     },
  { "widthx"        , VAR(widthx)          , PARAM_DOUBLE(1.  * 1e-6)     },
  { "widthy"        , VAR(widthy)          , PARAM_DOUBLE(1.  * 1e-6)     },
  { "widthz"        , VAR(widthz)          , PARAM_DOUBLE(15. * 1e-9)     },
  { "rot"           , VAR(rot)             , PARAM_DOUBLE(0.)             },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(12.*1836.)      },
  {},
};

#undef VAR

static void
psc_case_thinfoil_create(struct psc_case *_case)
{
  struct psc_pulse_flattop prm = {
    .xm = .01   * 1e-6,
    .ym = .01   * 1e-6,
    .zm = -301. * 1e-6,
    .dxm = 4.   * 1e-6,
    .dym = 4.   * 1e-6,
    .dzm = .1   * 1e-6,
    .zb  = 300. * 1e-6,
  };
  psc.pulse_p_z1 = psc_pulse_flattop_create(&prm);
}

static void
psc_case_thinfoil_set_from_options(struct psc_case *_case)
{
  psc.prm.nmax = 10000;
  psc.prm.cpum = 25000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 1.0e20;
  psc.prm.n0 = 1.0e29;

  psc.prm.nicell = 1000;

  psc.domain.length[0] = 0.02 * 1e-6;
  psc.domain.length[1] = 0.02 * 1e-6;
  psc.domain.length[2] = 2.0  * 1e-6;

  psc.domain.gdims[0] = 1;
  psc.domain.gdims[1] = 1;
  psc.domain.gdims[2] = 500;

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
psc_case_thinfoil_init_field(struct psc_case *_case, mfields_base_t *flds)
{
  struct psc *psc = _case->psc;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dy = psc->dx[1], dz = psc->dx[2], dt = psc->dt;
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);
      
      // FIXME, why this time?
      F3_BASE(pf, EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, 0.*dt);
      F3_BASE(pf, HX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
    } foreach_3d_g_end;
  }
}

static void
psc_case_thinfoil_init_npt(struct psc_case *_case, int kind, double x[3], 
			    struct psc_particle_npt *npt)
{
  struct psc_case_thinfoil *thinfoil = mrc_to_subobj(_case, struct psc_case_thinfoil);

  real Te = thinfoil->Te, Ti = thinfoil->Ti;

  real ld = psc.coeff.ld;

  real x0 = thinfoil->x0 / ld;
  real y0 = thinfoil->y0 / ld;
  real z0 = thinfoil->z0 / ld;
  real Lx = thinfoil->Lx / ld;
  real Ly = thinfoil->Ly / ld;
  real Lz = thinfoil->Lz / ld;
  real widthx = thinfoil->widthx / ld;
  real widthy = thinfoil->widthy / ld;
  real widthz = thinfoil->widthz / ld;
  real rot = thinfoil->rot * 2.*M_PI / 360.;

  real xr = x[0];
  real yr = cos(rot) * (x[1]-y0) - sin(rot) * (x[2]-z0) + y0;
  real zr = cos(rot) * (x[2]-z0) + sin(rot) * (x[1]-y0) + z0;

  real argx = (fabs(xr-x0)-widthx)/Lx;
  real argy = (fabs(yr-y0)-widthy)/Ly;
  real argz = (fabs(zr-z0)-widthz)/Lz;
  if (argx > 200.0) argx = 200.0;
  if (argy > 200.0) argy = 200.0;
  if (argz > 200.0) argz = 200.0;

  real dens = 1. / ((1. + exp(argx)) * (1. + exp(argy)) * (1. + exp(argz)));

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
    npt->m = thinfoil->mass_ratio;
    npt->n = dens;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_thinfoil_ops = {
  .name             = "thinfoil",
  .size             = sizeof(struct psc_case_thinfoil),
  .param_descr      = psc_case_thinfoil_descr,
  .create           = psc_case_thinfoil_create,
  .set_from_options = psc_case_thinfoil_set_from_options,
  .init_field       = psc_case_thinfoil_init_field,
  .init_npt         = psc_case_thinfoil_init_npt,
};
