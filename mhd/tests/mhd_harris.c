
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_crds_private.h>
#include <ggcm_mhd_crds_gen.h>
#include <ggcm_mhd_bnd.h>
#include <ggcm_mhd_diag.h>

#include <mrc_fld_as_double.h>
#include <mrc_domain.h>

#include <math.h>


#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>


// ======================================================================
// ggcm_mhd_ic subclass "harris"

struct ggcm_mhd_ic_harris {
  float n_0; //< peek density - n_inf
  float n_inf; //< background density
  float B_0; //< B field strength
  float cs_width; // current sheet width
  float pert; // strength of \psi perturbation (flux function)
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_run

static void
ggcm_mhd_ic_harris_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_harris *ic_harris = mrc_to_subobj(ic, struct ggcm_mhd_ic_harris);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  float cs_width = ic_harris->cs_width;

  float xl[3], xh[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  float lx = xh[0] - xl[0];
  float ly = xh[1] - xl[1];
  float x0 = xl[0];
  float y0 = xl[1];
  float kx = 2.0 * M_PI / lx;
  float ky = M_PI / ly;
  // u_th = constant thermal energy density (p / n)
  float u_th = sqr(ic_harris->B_0) / (2.0 * ic_harris->n_0);

  struct mrc_fld *fld_psi = mrc_domain_fld_create(fld->_domain, SW_2, "psi");
  mrc_fld_set_type(fld_psi, FLD_TYPE);
  mrc_fld_setup(fld_psi);
  mrc_fld_foreach(fld, ix,iy,iz, 2, 2) {
    float x = MRC_CRDX(crds, ix);
    float y = MRC_CRDY(crds, iy);

    //    A[2] = lx / (4*pi) * (1 - cos(2*kx*X)) * cos(ky*Y)
    // F3(fld_psi, 0, ix,iy,iz) = ic_harris->pert * ly / (4. * M_PI) * (1. - cos(2*ky*(y - y0))) * cos(kx*(x - x0));
    // taken from Birn et. al. 2001
    F3(fld_psi, 0, ix,iy,iz) = ic_harris->pert * \
                                   cos(ky*(y - y0 - ly/2.0)) * cos(kx*(x - x0 - lx/2.0));
  } mrc_fld_foreach_end;

  mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
    float y = MRC_CRDY(crds, iy);
    // eq
    RR(fld, ix,iy,iz) = ic_harris->n_inf + \
                         ic_harris->n_0 / (sqr(cosh((y - y0 - ly/2.0) / cs_width)));
    PP(fld, ix,iy,iz) = u_th * RR(fld, ix,iy,iz);
    BX(fld, ix,iy,iz) = ic_harris->B_0 * tanh((y - y0 - ly/2.0) / cs_width);
  } mrc_fld_foreach_end;

  // FIXME for nonuniform
  mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
    // perturbation from psi
    BX(fld, ix,iy,iz) +=
      (F3(fld_psi,0, ix,iy+1,iz) - F3(fld_psi,0, ix,iy,iz)) /
      (MRC_CRDY(crds,iy+1) - MRC_CRDY(crds, iy));
    BY(fld, ix,iy,iz) = -
      (F3(fld_psi,0, ix+1,iy,iz) - F3(fld_psi,0, ix,iy,iz)) /
      (MRC_CRDX(crds,ix+1) - MRC_CRDX(crds, ix));
  } mrc_fld_foreach_end;

  mrc_fld_destroy(fld_psi);
  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_harris, x)
static struct param ggcm_mhd_ic_harris_descr[] = {
  { "n_0"             , VAR(n_0)             , PARAM_FLOAT(1.0)         },
  { "n_inf"           , VAR(n_inf)           , PARAM_FLOAT(0.2)         },
  { "B_0"             , VAR(B_0)             , PARAM_FLOAT(1.0)         },
  { "cs_width"        , VAR(cs_width)        , PARAM_FLOAT(0.5)         },
  { "pert"            , VAR(pert)            , PARAM_FLOAT(0.1)         },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_harris_ops = {
  .name        = "harris",
  .size        = sizeof(struct ggcm_mhd_ic_harris),
  .param_descr = ggcm_mhd_ic_harris_descr,
  .run         = ggcm_mhd_ic_harris_run,
};

// ======================================================================
// ggcm_mhd_ic subclass "asymharris"

struct ggcm_mhd_ic_asymharris {
  float beta_min; // min plasma beta
  float n01; //< asymtotic density
  float n02; //< asymtotic density
  float B01; //< B field strength
  float B02; //< B field strength
  float cs_width; // current sheet width
  float pert; // strength of \psi perturbation (flux function)
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_asymharris_run

static void
ggcm_mhd_ic_asymharris_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_asymharris *ic_asymharris = mrc_to_subobj(ic, struct ggcm_mhd_ic_asymharris);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  float cs_width = ic_asymharris->cs_width;

  float xl[3], xh[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  float lx = xh[0] - xl[0];
  float ly = xh[1] - xl[1];
  float x0 = xl[0];
  float y0 = xl[1];
  float kx = 2.0 * M_PI / lx;
  float ky = M_PI / ly;
  float bmax = fmax(ic_asymharris->B01, ic_asymharris->B02);
  float c = (0.5 * sqr(bmax)) * (1.0 + ic_asymharris->beta_min);
  float n01 = ic_asymharris->n01;
  float n02 = ic_asymharris->n02;

  // mprintf("c = %f\n", c);

  struct mrc_fld *fld_psi = mrc_domain_fld_create(fld->_domain, SW_2, "psi");
  mrc_fld_setup(fld_psi);
  mrc_fld_foreach(fld, ix,iy,iz, 2, 2) {
    float x = MRC_CRDX(crds, ix);
    float y = MRC_CRDY(crds, iy);

    //    A[2] = lx / (4*pi) * (1 - cos(2*kx*X)) * cos(ky*Y)
    // F3(fld_psi, 0, ix,iy,iz) = ic_harris->pert * ly / (4. * M_PI) * (1. - cos(2*ky*(y - y0))) * cos(kx*(x - x0));
    // taken from Birn et. al. 2001
    F3(fld_psi, 0, ix,iy,iz) = ic_asymharris->pert * \
                                   cos(ky*(y - y0 - ly/2.0)) * cos(kx*(x - x0 - lx/2.0));
  } mrc_fld_foreach_end;

  mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
    float y = MRC_CRDY(crds, iy);
    float yprime = y - y0 - ly/2.0;  // y shifted to center of box
    float B0;

    if (yprime > 0.0) {
      B0 = ic_asymharris->B01;
    } else {
      B0 = ic_asymharris->B02;
    }

    BX(fld, ix,iy,iz) = B0 * tanh(yprime / cs_width);

    RR(fld, ix,iy,iz) = (0.5 * (n01 + n02) +
			 0.5 * (n01 - n02) * tanh(yprime / cs_width));
    PP(fld, ix,iy,iz) = c - (0.5 * sqr(BX(fld, ix,iy,iz)));

    // RR(fld, ix,iy,iz) = 0.2 +
    //                     n01 / (sqr(cosh(yprime / cs_width)));
    // PP(fld, ix,iy,iz) = u_th * RR(fld, ix,iy,iz);

  } mrc_fld_foreach_end;

  // FIXME for nonuniform
  mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
    // perturbation from psi
    BX(fld, ix,iy,iz) +=
      (F3(fld_psi,0, ix,iy+1,iz) - F3(fld_psi,0, ix,iy,iz)) /
      (MRC_CRDY(crds,iy+1) - MRC_CRDY(crds, iy));
    BY(fld, ix,iy,iz) = -
      (F3(fld_psi,0, ix+1,iy,iz) - F3(fld_psi,0, ix,iy,iz)) /
      (MRC_CRDX(crds,ix+1) - MRC_CRDX(crds, ix));
  } mrc_fld_foreach_end;

  mrc_fld_destroy(fld_psi);
  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_asymharris_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_asymharris, x)
static struct param ggcm_mhd_ic_asymharris_descr[] = {
  { "beta_min"        , VAR(beta_min)        , PARAM_FLOAT(4.0)         },
  { "n01"             , VAR(n01)             , PARAM_FLOAT(1.0)         },
  { "n02"             , VAR(n02)             , PARAM_FLOAT(1.0)         },
  { "B01"             , VAR(B01)             , PARAM_FLOAT(1.0)         },
  { "B02"             , VAR(B02)             , PARAM_FLOAT(1.0)         },
  { "cs_width"        , VAR(cs_width)        , PARAM_FLOAT(0.5)         },
  { "pert"            , VAR(pert)            , PARAM_FLOAT(0.1)         },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_asymharris_ops = {
  .name        = "asymharris",
  .size        = sizeof(struct ggcm_mhd_ic_asymharris),
  .param_descr = ggcm_mhd_ic_asymharris_descr,
  .run         = ggcm_mhd_ic_asymharris_run,
};


// ======================================================================
// ggcm_mhd subclass "ot"

// ----------------------------------------------------------------------
// ggcm_mhd_ot_create

static void
ggcm_mhd_harris_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  -12.8, -6.4, 0.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {   12.8,  6.4, 0.1 });
  ggcm_mhd_set_param_float(mhd, "diffconstant", 0.005);

  ggcm_mhd_bnd_set_type(mhd->bnd, "conducting_y_c2");
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcy", BC_NONE);
  mrc_domain_set_param_int(mhd->domain, "bcz", BC_PERIODIC);

  ggcm_mhd_ic_set_type(mhd->ic, "harris");
  ggcm_mhd_step_set_type(mhd->step, "c2_float");
}

// ----------------------------------------------------------------------
// ggcm_mhd_ot_ops

static struct ggcm_mhd_ops ggcm_mhd_harris_ops = {
  .name             = "harris",
  .create           = ggcm_mhd_harris_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_harris_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_harris_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_asymharris_ops);

  return ggcm_mhd_main(&argc, &argv);
}
