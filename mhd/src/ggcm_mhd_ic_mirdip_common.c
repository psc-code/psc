#include <stdlib.h>
#include <math.h>

#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bndsw.h"
#include "ggcm_mhd_dipole.h"
#include "ggcm_mhd_crds.h"

#include "mrc_domain.h"

// ======================================================================
// ggcm_mhd_ic subclass "mirdip"

struct ggcm_mhd_ic_mirdip {
  float xxx1;
  float xxx2;
  float xmir;
  float rrini;
  float prat; // determines initial interior pressure as a fraction of solar wind pressure
  float stretch_tail; // vx == 0 is stretched by this factor on the tail side

  float dipole_moment[3];

  // solar wind values to use if there is not "bndsw" object around
  double bnvals[SW_NR];
};

#define ggcm_mhd_ic_mirdip(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_mirdip)

// ----------------------------------------------------------------------
// get_solar_wind

static void
get_solar_wind(struct ggcm_mhd_ic *ic, float vals[])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  struct ggcm_mhd *mhd = ic->mhd;

  vals[SW_RR] = sub->bnvals[SW_RR] / mhd->rrnorm;
  vals[SW_VX] = sub->bnvals[SW_VX] / mhd->vvnorm;
  vals[SW_VY] = sub->bnvals[SW_VY] / mhd->vvnorm;
  vals[SW_VZ] = sub->bnvals[SW_VZ] / mhd->vvnorm;
  vals[SW_PP] = sub->bnvals[SW_PP] / mhd->ppnorm;
  vals[SW_BX] = sub->bnvals[SW_BX] / mhd->bbnorm;
  vals[SW_BY] = sub->bnvals[SW_BY] / mhd->bbnorm;
  vals[SW_BZ] = sub->bnvals[SW_BZ] / mhd->bbnorm;
}

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
  if (xmir != 0. && x < xmir) {
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
  mrc_fld_data_t rrini = sub->rrini, stretch_tail = sub->stretch_tail;

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 2, 2) {
      mrc_fld_data_t tmplam = lmbda(MRC_MCRDX(crds, ix, p), -xxx1, -xxx2);
      RR_(fld, ix,iy,iz, p) = 
	tmplam * vals[SW_RR] + (1. - tmplam) * (rrini / mhd->rrnorm);
      
      VX_(fld, ix,iy,iz, p) = vxsta1(MRC_MCRDX(crds, ix, p), MRC_MCRDY(crds, iy, p), MRC_MCRDZ(crds, iz, p),
				     vals[SW_VX], xxx2, xxx1, stretch_tail, xmir);
      VY_(fld, ix,iy,iz, p) = 0.;
      VZ_(fld, ix,iy,iz, p) = 0.;
      
      mrc_fld_data_t ppmin = sub->prat * vals[SW_PP];
      PP_(fld, ix,iy,iz, p) = vals[SW_PP] * tmplam + (1. - tmplam) * ppmin;
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_get_mhd_dipole

static struct ggcm_mhd_dipole *
ggcm_mhd_ic_mirdip_get_mhd_dipole(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_dipole *mhd_dipole = ggcm_mhd_get_var_obj(ic->mhd, "mhd_dipole");
  if (mhd_dipole) {
    ggcm_mhd_dipole_get(mhd_dipole);
  } else {
    mhd_dipole = ggcm_mhd_dipole_create(ggcm_mhd_ic_comm(ic));
    ggcm_mhd_dipole_set_type(mhd_dipole, FLD_TYPE);
    ggcm_mhd_dipole_set_from_options(mhd_dipole);
    ggcm_mhd_dipole_set_param_obj(mhd_dipole, "mhd", ic->mhd);
    ggcm_mhd_dipole_setup(mhd_dipole);
  }

  return mhd_dipole;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_vector_potential_bg

static double
ggcm_mhd_ic_mirdip_vector_potential_bg(struct ggcm_mhd_ic *ic, int m, double x[3])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  static struct ggcm_mhd_dipole *mhd_dipole;
  static float vals[SW_NR];
  static bool first_time = true;
  if (first_time) {
    mhd_dipole = ggcm_mhd_ic_mirdip_get_mhd_dipole(ic);
    get_solar_wind(ic, vals);
    first_time = false;
  }

  // get main dipole vector potential
  double A = ggcm_mhd_dipole_vector_potential(mhd_dipole, m, x, (float [3]) { 0.f, 0.f, 0.f },
					      sub->dipole_moment, sub->xmir);

  // add IMF vector potential
  switch (m) {
  case 1: A += vals[SW_BZ] * x[0]; break;
  case 2: A += vals[SW_BX] * x[1] - vals[SW_BY] * x[0]; break;
  }

  return A;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_vector_potential

static double
ggcm_mhd_ic_mirdip_vector_potential(struct ggcm_mhd_ic *ic, int m, double x[3])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  static struct ggcm_mhd_dipole *mhd_dipole;
  static float vals[SW_NR];
  static bool first_time = true;
  if (first_time) {
    mhd_dipole = ggcm_mhd_ic_mirdip_get_mhd_dipole(ic);
    get_solar_wind(ic, vals);
    first_time = false;
  }

  // get main dipole vector potential
  double A = ggcm_mhd_dipole_vector_potential(mhd_dipole, m, x, (float [3]) { 0.f, 0.f, 0.f },
					      sub->dipole_moment, sub->xmir);

  if (sub->xmir != 0.0) {
    // add mirror dipole vector potential
    sub->dipole_moment[0] *= -1.f;
    A += ggcm_mhd_dipole_vector_potential(mhd_dipole, m, x, (float [3]) { 2.f * sub->xmir, 0.f, 0.f },
					  sub->dipole_moment, sub->xmir);
    sub->dipole_moment[0] *= -1.f;
  }

  // add IMF vector potential
  switch (m) {
  case 1: A += vals[SW_BZ] * x[0]; break;
  case 2: A += vals[SW_BX] * x[1] - vals[SW_BY] * x[0]; break;
  }

  // subtract out B0-vector potential
  // this is somewhat odd, as in, why are we adding the main dipole and then subtract it here again,
  // but the main dipole field is actually not cut-off at xmir, whereas this one here is, so they
  // different
  A -= ggcm_mhd_ic_mirdip_vector_potential_bg(ic, m, x);

  return A;
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
  int mhd_type_save;
  mrc_fld_get_param_int(ic->mhd->fld, "mhd_type", &mhd_type_save);
  int mhd_type = MT_PRIMITIVE;
  if (mhd_type_save == MT_FULLY_CONSERVATIVE_CC) {
    mhd_type = MT_PRIMITIVE_CC;
  }
  mrc_fld_set_param_int(ic->mhd->fld, "mhd_type", mhd_type);

  float vals[SW_NR];
  get_solar_wind(ic, vals);

  ggcm_mhd_ic_mirdip_ini1(ic, vals);

  mrc_fld_set_param_int(ic->mhd->fld, "mhd_type", mhd_type_save);
  if (mhd_type == MT_PRIMITIVE) {
    ggcm_mhd_convert_from_primitive(ic->mhd, ic->mhd->fld);
  } else if (mhd_type == MT_PRIMITIVE_CC) {
    ggcm_mhd_convert_from_primitive_cc(ic->mhd, ic->mhd->fld);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_mirdip, x)
static struct param ggcm_mhd_ic_mirdip_descr[] = {
  { "xxx1"         , VAR(xxx1)         , PARAM_FLOAT(14.0)         },
  { "xxx2"         , VAR(xxx2)         , PARAM_FLOAT(12.0)         },
  { "xmir"         , VAR(xmir)         , PARAM_FLOAT(-15.0)        },  // off if == 0.0
  { "rrini"        , VAR(rrini)        , PARAM_FLOAT(3.0)          },
  { "prat"         , VAR(prat)         , PARAM_FLOAT(.5)           },
  { "stretch_tail" , VAR(stretch_tail) , PARAM_FLOAT(2.)           },

  { "dipole_moment", VAR(dipole_moment), PARAM_FLOAT3(0., 0., -1.) },

  { "rr"           , VAR(bnvals[SW_RR]), PARAM_DOUBLE(1.)          },
  { "pp"           , VAR(bnvals[SW_PP]), PARAM_DOUBLE(1.)          },
  { "vx"           , VAR(bnvals[SW_VX]), PARAM_DOUBLE(0.)          },
  { "vy"           , VAR(bnvals[SW_VY]), PARAM_DOUBLE(0.)          },
  { "vz"           , VAR(bnvals[SW_VZ]), PARAM_DOUBLE(0.)          },
  { "bx"           , VAR(bnvals[SW_BX]), PARAM_DOUBLE(0.)          },
  { "by"           , VAR(bnvals[SW_BY]), PARAM_DOUBLE(0.)          },
  { "bz"           , VAR(bnvals[SW_BZ]), PARAM_DOUBLE(0.)          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic subclass "mirdip"

struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip_ops = {
  .name                = ggcm_mhd_ic_mirdip_name,
  .size                = sizeof(struct ggcm_mhd_ic_mirdip),
  .param_descr         = ggcm_mhd_ic_mirdip_descr,
  .run                 = ggcm_mhd_ic_mirdip_run,
  .vector_potential_bg = ggcm_mhd_ic_mirdip_vector_potential_bg,
  .vector_potential    = ggcm_mhd_ic_mirdip_vector_potential,
};
