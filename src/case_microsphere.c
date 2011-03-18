
#include "psc.h"
#include "psc_case_private.h"
#include "psc_push_fields.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ----------------------------------------------------------------------
// case microsphere
//

struct psc_case_microsphere {
  // parameters
};

#define VAR(x) (void *)offsetof(struct psc_case_microsphere, x)

static struct param psc_case_microsphere_descr[] = {
  {},
};

#undef VAR

static double
microsphere_dens(struct psc_case *_case, double x[3])
{
  // center
  const double xc[3] = { 5. * 1e-6, 5. * 1e-6, 5. * 1e-6 };
  const double radius = 2. * 1e-6;
  const double thickness = .4 * 1e-6;

  double xr[3];
  for (int d = 0; d < 3; d++) {
    xr[d] = x[d] * psc.coeff.ld - xc[0];
  };

  double r = sqrt(sqr(xr[0]) + sqr(xr[1]) + sqr(xr[2]));

  return exp(-sqr((r - radius) / thickness));
}

static void
psc_case_microsphere_set_from_options(struct psc_case *_case)
{
  //  struct psc_case_microsphere *msphere = mrc_to_subobj(_case, struct psc_case_microsphere);

  psc.prm.nicell = 10;
  psc.prm.nmax = 201;

  psc.domain.length[0] = 10 * 1e-6;
  psc.domain.length[1] = 10 * 1e-6;
  psc.domain.length[2] = 10 * 1e-6;

  psc.domain.gdims[0] = 128;
  psc.domain.gdims[1] = 128;
  psc.domain.gdims[2] = 128;

  psc.domain.bnd_fld_lo[0] = BND_FLD_OPEN;
  psc.domain.bnd_fld_hi[0] = BND_FLD_OPEN;
  psc.domain.bnd_fld_lo[1] = BND_FLD_OPEN;
  psc.domain.bnd_fld_hi[1] = BND_FLD_OPEN;
  psc.domain.bnd_fld_lo[2] = BND_FLD_OPEN;
  psc.domain.bnd_fld_hi[2] = BND_FLD_OPEN;
  psc.domain.bnd_part[0] = BND_PART_PERIODIC; // FIXME
  psc.domain.bnd_part[1] = BND_PART_PERIODIC;
  psc.domain.bnd_part[2] = BND_PART_PERIODIC;

  psc_push_fields_set_type(psc.push_fields, "fortran");

  double *length = psc.domain.length;
  double width_normal = 2. * 1e-6;
  double width_par    = 3. * 1e-6;

  struct psc_pulse_gauss prm_p_x1 = {
    .xm  = 0. * length[0],
    .ym  = .5 * length[1],
    .zm  = .5 * length[2],
    .dxm = width_normal,
    .dym = width_par,
    .dzm = width_par,
    .k  = { 1., 0., 0. },
  };

  struct psc_pulse_gauss prm_p_x2 = {
    .xm  = 1. * length[0],
    .ym  = .5 * length[1],
    .zm  = .5 * length[2],
    .dxm = width_normal,
    .dym = width_par,
    .dzm = width_par,
    .k  = { -1., 0., 0. },
  };

  psc.pulse_p_x1 = psc_pulse_gauss_create(&prm_p_x1);
  psc.pulse_p_x2 = psc_pulse_gauss_create(&prm_p_x2);

  struct psc_pulse_gauss prm_p_y1 = {
    .xm  = .5 * length[0],
    .ym  = 0. * length[1],
    .zm  = .5 * length[2],
    .dxm = width_par,
    .dym = width_normal,
    .dzm = width_par,
    .k  = { 0., 1., 0. },
  };

  struct psc_pulse_gauss prm_p_y2 = {
    .xm  = .5 * length[0],
    .ym  = 1. * length[1],
    .zm  = .5 * length[2],
    .dxm = width_par,
    .dym = width_normal,
    .dzm = width_par,
    .k  = { 0., -1., 0. },
  };

  psc.pulse_p_y1 = psc_pulse_gauss_create(&prm_p_y1);
  psc.pulse_p_y2 = psc_pulse_gauss_create(&prm_p_y2);

  struct psc_pulse_gauss prm_p_z1 = {
    .xm  = .5 * length[0],
    .ym  = .5 * length[1],
    .zm  = 0. * length[2],
    .dxm = width_par,
    .dym = width_par,
    .dzm = width_normal,
    .k  = { 0., 0., 1. },
  };

  struct psc_pulse_gauss prm_p_z2 = {
    .xm  = .5 * length[0],
    .ym  = .5 * length[1],
    .zm  = 1. * length[2],
    .dxm = width_par,
    .dym = width_par,
    .dzm = width_normal,
    .k  = { 0., 0., -1. },
  };

  psc.pulse_p_z1 = psc_pulse_gauss_create(&prm_p_z1);
  psc.pulse_p_z2 = psc_pulse_gauss_create(&prm_p_z2);
}

static void
psc_case_microsphere_init_field(struct psc_case *_case, mfields_base_t *flds)
{
  //  struct psc_case_microsphere *msphere = mrc_to_subobj(_case, struct psc_case_microsphere);
  struct psc *psc = _case->psc;

  // FIXME, do we need the ghost points?
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    psc_foreach_3d_g(psc, p, jx, jy, jz) {
      double dx = psc->dx[0], dy = psc->dx[1], dz = psc->dx[2];
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);

      if (psc->pulse_p_x1) {
	F3_BASE(pf, EZ, jx,jy,jz) +=  psc_p_pulse_x1(xx        , yy, zz + .5*dz, 0.);
	F3_BASE(pf, HY, jx,jy,jz) += -psc_p_pulse_x1(xx + .5*dx, yy, zz + .5*dz, 0.);
      }
      
      if (psc->pulse_p_x2) {
	F3_BASE(pf, EZ, jx,jy,jz) +=  psc_p_pulse_x2(xx        , yy, zz + .5*dz, 0.);
	F3_BASE(pf, HY, jx,jy,jz) +=  psc_p_pulse_x2(xx + .5*dx, yy, zz + .5*dz, 0.);
      }
      
      if (psc->pulse_p_y1) {
	F3_BASE(pf, EZ, jx,jy,jz) +=  psc_p_pulse_y1(xx, yy        , zz + .5*dz, 0.);
	F3_BASE(pf, HX, jx,jy,jz) +=  psc_p_pulse_y1(xx, yy + .5*dy, zz + .5*dz, 0.);
      }
      
      if (psc->pulse_p_y2) {
	F3_BASE(pf, EZ, jx,jy,jz) +=  psc_p_pulse_y2(xx, yy        , zz + .5*dz, 0.);
	F3_BASE(pf, HX, jx,jy,jz) += -psc_p_pulse_y2(xx, yy + .5*dy, zz + .5*dz, 0.);
      }
      
      if (psc->pulse_p_z1) {
	F3_BASE(pf, EY, jx,jy,jz) +=  psc_p_pulse_z1(xx, yy + .5*dy, zz        , 0.);
	F3_BASE(pf, HX, jx,jy,jz) += -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.);
      }
      
      if (psc->pulse_p_z2) {
	F3_BASE(pf, EY, jx,jy,jz) +=  psc_p_pulse_z2(xx, yy + .5*dy, zz        , 0.);
	F3_BASE(pf, HX, jx,jy,jz) +=  psc_p_pulse_z2(xx, yy + .5*dy, zz + .5*dz, 0.);
      }
    } psc_foreach_3d_g_end;
  }
}

static void
psc_case_microsphere_init_npt(struct psc_case *_case, int kind, double x[3],
			     struct psc_particle_npt *npt)
{
  //  struct psc_case_microsphere *msphere = mrc_to_subobj(_case, struct psc_case_microsphere);
  double dens = microsphere_dens(_case, x);

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1;
    npt->n = dens;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 100;
    npt->n = dens;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_microsphere_ops = {
  .name             = "microsphere",
  .size             = sizeof(struct psc_case_microsphere),
  .param_descr      = psc_case_microsphere_descr,
  .set_from_options = psc_case_microsphere_set_from_options,
  .init_field       = psc_case_microsphere_init_field,
  .init_npt         = psc_case_microsphere_init_npt,
};

