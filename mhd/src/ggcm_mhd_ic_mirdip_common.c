#include <stdlib.h>
#include <math.h>

#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bndsw.h"
#include "ggcm_mhd_dipole.h"

#include "mrc_domain.h"

// ======================================================================
// ggcm_mhd_ic subclass "mirdip"

struct ggcm_mhd_ic_mirdip {
  float xxx1;
  float xxx2;
  float xmir;
  float rrini;
  float r1lim;
  float prat; // determines initial interior pressure as a fraction of solar wind pressure

  float dipole_moment[3];
};

#define ggcm_mhd_ic_mirdip(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_mirdip)

// ----------------------------------------------------------------------
// lmbda

static inline mrc_fld_data_t
lmbda(mrc_fld_data_t XX, mrc_fld_data_t R1, mrc_fld_data_t R2)
{
  mrc_fld_data_t LL = (XX - R2) / (R1 - R2);
  LL = fmax(LL, 0.);
  LL = fmin(LL, 1.);
  return LL;
}

// ----------------------------------------------------------------------
// vxsta1

static inline mrc_fld_data_t
vxsta1(mrc_fld_data_t x, mrc_fld_data_t y, mrc_fld_data_t z, mrc_fld_data_t v0,
       mrc_fld_data_t r1, mrc_fld_data_t r2, mrc_fld_data_t str, mrc_fld_data_t xmir)
{
  mrc_fld_data_t xx = x;
  if (xx > 0.) {
    xx = x / str;
  }
  mrc_fld_data_t r = sqrt(sqr(xx) + sqr(y) + sqr(z));
  mrc_fld_data_t s = (r - r1) / (r2 - r1);
  s = fmax(s, 0.);
  s = fmin(s, 1.);
  if (x < xmir) {
    return v0;
  } else {
    return v0 * s;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_ini1

static void
ggcm_mhd_ic_mirdip_ini1(struct ggcm_mhd_ic *ic, float vals[])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);

  mrc_fld_data_t xxx1 = sub->xxx1, xxx2 = sub->xxx2, xmir = sub->xmir;
  mrc_fld_data_t rrini = sub->rrini;

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 2, 2) {
      mrc_fld_data_t tmplam = lmbda(MRC_MCRDX(crds, ix, p), -xxx1, -xxx2);
      RR_(fld, ix,iy,iz, p) = 
	tmplam * vals[SW_RR] + (1. - tmplam) * (rrini / mhd->par.rrnorm);
      
      VX_(fld, ix,iy,iz, p) = vxsta1(MRC_MCRDX(crds, ix, p), MRC_MCRDY(crds, iy, p), MRC_MCRDZ(crds, iz, p),
				     vals[SW_VX], xxx2, xxx1, 2., xmir);
      VY_(fld, ix,iy,iz, p) = 0.;
      VZ_(fld, ix,iy,iz, p) = 0.;
      
      mrc_fld_data_t ppmin = sub->prat * vals[SW_PP];
      PP_(fld, ix,iy,iz, p) = vals[SW_PP] * tmplam + (1. - tmplam) * ppmin;
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_ini_b

static void
ggcm_mhd_ic_mirdip_ini_b(struct ggcm_mhd_ic *ic, float b_sw[3])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  struct ggcm_mhd_dipole *mhd_dipole = ggcm_mhd_get_var_obj(ic->mhd, "mhd_dipole");
  if (!mhd_dipole) {
    return;
  }

  float x0[3] = {0.0, 0.0, 0.0};

  /* mprintf("mirdip_ini_b dipole moment: %f %f %f\n", sub->dipole_moment[0], */
  /* 	  sub->dipole_moment[1], sub->dipole_moment[2]); */

  struct mrc_fld *_b = mrc_fld_make_view(ic->mhd->fld, BX, BX + 3);
  ggcm_mhd_dipole_add_dipole(mhd_dipole, _b, x0, sub->dipole_moment, sub->xmir, 0.);

  if (sub->xmir != 0.0) {
    x0[0] = 2.0 * sub->xmir;
    sub->dipole_moment[0] *= -1.0;
    ggcm_mhd_dipole_add_dipole(mhd_dipole, _b, x0, sub->dipole_moment, sub->xmir, 1.);
    sub->dipole_moment[0] *= -1.0;
  }

  struct mrc_fld *b = mrc_fld_get_as(_b, FLD_TYPE);

  // finish up
  for (int p = 0; p < mrc_fld_nr_patches(b); p++) {
    mrc_fld_foreach(b, ix,iy,iz, 0, 1) {
      for (int i = 0; i < 3; i++){
	float crd_fc[3];
	ggcm_mhd_get_crds_fc(ic->mhd, ix,iy,iz, p, i, crd_fc);
	
	// zero B if inside r1lim
	float r = sqrtf(sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]));
	if (r < sub->r1lim) {
	  M3(b, i, ix,iy,iz, p) = 0.0;
	}
	
	// add B_IMF
	M3(b, i, ix,iy,iz, p) += b_sw[i];
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(b, _b);

  mrc_fld_destroy(_b);

}

// ----------------------------------------------------------------------
// ggcm_mhd_mirdip_ic_run

static void
ggcm_mhd_ic_mirdip_run(struct ggcm_mhd_ic *ic)
{
  // ini_b is done when we still have the primitive variables stored in
  // mhd->fld, so it needs to use the right ("regular") staggering because
  // the convert_from_primitive() will fix up the staggering to whatever
  // we really need in the end.
  int mhd_type;
  mrc_fld_get_param_int(ic->mhd->fld, "mhd_type", &mhd_type);
  mrc_fld_set_param_int(ic->mhd->fld, "mhd_type", MT_PRIMITIVE);

  float vals[SW_NR];
  struct ggcm_mhd_bndsw *bndsw = ggcm_mhd_get_var_obj(ic->mhd, "bndsw");
  ggcm_mhd_bndsw_get_initial(bndsw, vals);
  ggcm_mhd_ic_mirdip_ini1(ic, vals);
  ggcm_mhd_ic_mirdip_ini_b(ic, &vals[SW_BX]);

  mrc_fld_set_param_int(ic->mhd->fld, "mhd_type", mhd_type);
  ggcm_mhd_convert_from_primitive(ic->mhd, ic->mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_mirdip_ic_init_ymask

static void
ggcm_mhd_ic_mirdip_init_ymask(struct ggcm_mhd_ic *ic, struct mrc_fld *ymask)
{
  ymaskn_c(ic->mhd, ymask);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_mirdip, x)
static struct param ggcm_mhd_ic_mirdip_descr[] = {
  { "xxx1"         , VAR(xxx1)         , PARAM_FLOAT(14.0)         },
  { "xxx2"         , VAR(xxx2)         , PARAM_FLOAT(12.0)         },
  { "xmir"         , VAR(xmir)         , PARAM_FLOAT(-15.0)        },  // off if == 0.0
  { "rrini"        , VAR(rrini)        , PARAM_FLOAT(3.0)          },
  { "r1lim"        , VAR(r1lim)        , PARAM_FLOAT(1.5)          },
  { "prat"         , VAR(prat)         , PARAM_FLOAT(.5)           },

  { "dipole_moment", VAR(dipole_moment), PARAM_FLOAT3(0., 0., -1.) },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic subclass "mirdip"

struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip_ops = {
  .name             = ggcm_mhd_ic_mirdip_name,
  .size             = sizeof(struct ggcm_mhd_ic_mirdip),
  .param_descr      = ggcm_mhd_ic_mirdip_descr,
  .run              = ggcm_mhd_ic_mirdip_run,
  .init_ymask       = ggcm_mhd_ic_mirdip_init_ymask,
};
