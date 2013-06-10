
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_flds.h"
#include "ggcmtest.h"

#include <mrc_domain.h>
#include <mrc_crds.h>
#include <math.h>

// ======================================================================
// ggcm_mhd_ic subclass "cos"

struct ggcm_mhd_ic_cos {
  float n_b; //< background density

  float pert;
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_cos_run

static void
ggcm_mhd_ic_cos_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_cos *ic_cos = mrc_to_subobj(ic, struct ggcm_mhd_ic_cos);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, "fortran");
  struct mrc_crds *crds = mrc_domain_get_crds(f->domain);
  float n_b = ic_cos->n_b;

  float xl[3], xh[3];
  mrc_crds_get_xl_xh(crds, xl, xh);
  real lx = xh[0] - xl[0], ly = xh[1] - xl[1];
  real kx = 2. * M_PI / lx, ky = 2. * M_PI / ly;

  struct mrc_fld *fld_psi = mrc_domain_fld_create(f->domain, SW_2);
  mrc_fld_setup(fld_psi);
  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    real x = MRC_CRDX(crds, ix), y = MRC_CRDY(crds, iy);

    //    A[2] = lx / (4*pi) * (1 - cos(2*kx*X)) * sin(ky*Y)
    MRC_F3(fld_psi, 0, ix,iy,iz) = ic_cos->pert * lx / (4. * M_PI) * (1. - cos(2*kx*x)) * sin(ky*y);
  } mrc_fld_foreach_end;

  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    real x = MRC_CRDX(crds, ix);
    // eq
    MRC_F3(f, _RR, ix,iy,iz) = n_b + .5 * sqr(sin(kx * x));
    MRC_F3(f, _PP, ix,iy,iz) = n_b + .5 * sqr(sin(kx * x));
    MRC_F3(f, _BY, ix,iy,iz) = cos(kx * x);
  } mrc_fld_foreach_end;

  ggcmtest_init_from_primitive(mhd);

  // FIXME for nonuniform
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    // perturbation from psi
    MRC_F3(f, _B1X, ix,iy,iz) +=
      (MRC_F3(fld_psi,0, ix+1,iy+1,iz) - MRC_F3(fld_psi,0, ix+1,iy,iz)) /
      (MRC_CRDY(crds,iy+1) - MRC_CRDY(crds, iy));
    MRC_F3(f, _B1Y, ix,iy,iz) += -
      (MRC_F3(fld_psi,0, ix+1,iy+1,iz) - MRC_F3(fld_psi,0, ix,iy+1,iz)) /
      (MRC_CRDX(crds,ix+1) - MRC_CRDX(crds, ix));
  } mrc_fld_foreach_end;

  mrc_fld_destroy(fld_psi);
  mrc_fld_put_as(f, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_cos_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_cos, x)
static struct param ggcm_mhd_ic_cos_descr[] = {
  { "n_b"             , VAR(n_b)             , PARAM_FLOAT(.2)         },
  { "pert"            , VAR(pert)            , PARAM_FLOAT(1e-2)       },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_cos_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_cos_ops = {
  .name        = "cos",
  .size        = sizeof(struct ggcm_mhd_ic_cos),
  .param_descr = ggcm_mhd_ic_cos_descr,
  .run         = ggcm_mhd_ic_cos_run,
};

