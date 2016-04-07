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

  // solar wind values to use if there is no "bndsw" object around
  double bnvals[SW_NR];

  // state
  double bnvals_code[SW_NR]; // normalized to code units
  double rrini_code;

  struct ggcm_mhd_dipole *mhd_dipole;
};

#define ggcm_mhd_ic_mirdip(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_mirdip)

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_setup

static void
ggcm_mhd_ic_mirdip_setup(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  struct ggcm_mhd *mhd = ic->mhd;

  sub->bnvals_code[SW_RR] = sub->bnvals[SW_RR] / mhd->rrnorm;
  sub->bnvals_code[SW_VX] = sub->bnvals[SW_VX] / mhd->vvnorm;
  sub->bnvals_code[SW_VY] = sub->bnvals[SW_VY] / mhd->vvnorm;
  sub->bnvals_code[SW_VZ] = sub->bnvals[SW_VZ] / mhd->vvnorm;
  sub->bnvals_code[SW_PP] = sub->bnvals[SW_PP] / mhd->ppnorm;
  sub->bnvals_code[SW_BX] = sub->bnvals[SW_BX] / mhd->bbnorm;
  sub->bnvals_code[SW_BY] = sub->bnvals[SW_BY] / mhd->bbnorm;
  sub->bnvals_code[SW_BZ] = sub->bnvals[SW_BZ] / mhd->bbnorm;

  sub->rrini_code = sub->rrini / mhd->rrnorm;

  // need a mhd_dipole object that'll take care of setting the dipole(s) involved
  sub->mhd_dipole = ggcm_mhd_get_var_obj(ic->mhd, "mhd_dipole");
  if (sub->mhd_dipole) {
    // if run as openggcm, ggcm_mhd has a mhd_dipole member, so use that
    ggcm_mhd_dipole_get(sub->mhd_dipole);
  } else {
    // otherwise, we gotta make our own
    sub->mhd_dipole = ggcm_mhd_dipole_create(ggcm_mhd_ic_comm(ic));
    ggcm_mhd_dipole_set_type(sub->mhd_dipole, FLD_TYPE);
    ggcm_mhd_dipole_set_from_options(sub->mhd_dipole);
    ggcm_mhd_dipole_set_param_obj(sub->mhd_dipole, "mhd", ic->mhd);
    ggcm_mhd_dipole_setup(sub->mhd_dipole);
  }
}

// ----------------------------------------------------------------------
// lmbda

static inline mrc_fld_data_t
lmbda(mrc_fld_data_t XX, mrc_fld_data_t R1, mrc_fld_data_t R2)
{
  mrc_fld_data_t LL = (XX - R2) / (R1 - R2);
  return fmin(fmax(LL, 0.f), 1.f);
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
  s = fmin(fmax(s, 0.), 1.f);
  if (xmir != 0. && x < xmir) {
    return v0;
  } else {
    return v0 * s;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_primitive

static double
ggcm_mhd_ic_mirdip_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);

  double xx = crd[0], yy = crd[1], zz = crd[2];

  double *vals = sub->bnvals_code;
  mrc_fld_data_t tmplam = lmbda(xx, -sub->xxx1, -sub->xxx2);

  switch (m) {
  case RR: return tmplam * vals[SW_RR] + (1.f - tmplam) * sub->rrini_code;
  case PP: return tmplam * vals[SW_PP] + (1.f - tmplam) * sub->prat * vals[SW_PP];
  case VX: return vxsta1(xx, yy, zz, vals[SW_VX], sub->xxx2, sub->xxx1, sub->stretch_tail, sub->xmir);
  case VY: return 0.;
  case VZ: return 0.;
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_vector_potential_bg

static double
ggcm_mhd_ic_mirdip_vector_potential_bg(struct ggcm_mhd_ic *ic, int m, double x[3])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);

  // get main dipole vector potential
  double A = ggcm_mhd_dipole_vector_potential(sub->mhd_dipole, m, x, (float [3]) { 0.f, 0.f, 0.f },
					      sub->dipole_moment, sub->xmir);

  // add IMF vector potential
  double *vals = sub->bnvals_code;
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

  // get main dipole vector potential
  double A = ggcm_mhd_dipole_vector_potential(sub->mhd_dipole, m, x, (float [3]) { 0.f, 0.f, 0.f },
					      sub->dipole_moment, sub->xmir);

  if (sub->xmir != 0.0) {
    // add mirror dipole vector potential
    sub->dipole_moment[0] *= -1.f;
    A += ggcm_mhd_dipole_vector_potential(sub->mhd_dipole, m, x, (float [3]) { 2.f * sub->xmir, 0.f, 0.f },
					  sub->dipole_moment, sub->xmir);
    sub->dipole_moment[0] *= -1.f;
  }

  // add IMF vector potential
  double *vals = sub->bnvals_code;
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
  .setup               = ggcm_mhd_ic_mirdip_setup,
  .primitive           = ggcm_mhd_ic_mirdip_primitive,
  .vector_potential    = ggcm_mhd_ic_mirdip_vector_potential,
  .vector_potential_bg = ggcm_mhd_ic_mirdip_vector_potential_bg,
};
