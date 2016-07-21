#include <stdlib.h>
#include <math.h>

#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bndsw.h"

#include "mrc_domain.h"

// ======================================================================
// ggcm_mhd_ic subclass "obstacle"

struct ggcm_mhd_ic_obstacle {
  float r1; // radius where velocity is ramped down all the way to 0
  float r2; // radius where velocity ramp down starts

  // solar wind values
  double bnvals[SW_NR];
};

#define ggcm_mhd_ic_obstacle(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_obstacle)

// ----------------------------------------------------------------------
// get_solar_wind

static void
get_solar_wind(struct ggcm_mhd_ic *ic, mrc_fld_data_t vals[])
{
  struct ggcm_mhd_ic_obstacle *sub = ggcm_mhd_ic_obstacle(ic);
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
// vxsta1

static inline mrc_fld_data_t
vxsta1(mrc_fld_data_t x, mrc_fld_data_t y, mrc_fld_data_t z, mrc_fld_data_t v0,
       mrc_fld_data_t r1, mrc_fld_data_t r2)
{
  mrc_fld_data_t r = mrc_fld_sqrt(sqr(x) + sqr(y) + sqr(z));
  mrc_fld_data_t s = (r - r1) / (r2 - r1);
  s = mrc_fld_min(mrc_fld_max(s, 0.), 1.);
  return v0 * s;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_obstacle_primitive

static double
ggcm_mhd_ic_obstacle_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_obstacle *sub = ggcm_mhd_ic_obstacle(ic);

  mrc_fld_data_t vals[SW_NR];
  get_solar_wind(ic, vals);

  double xx = crd[0], yy = crd[1], zz = crd[2];

  switch (m) {
  case RR: return vals[SW_RR];
  case PP: return vals[SW_PP];
  case VX: return vxsta1(xx, yy, zz, vals[SW_VX], sub->r1, sub->r2);
  case VY: return vxsta1(xx, yy, zz, vals[SW_VY], sub->r1, sub->r2);
  case VZ: return vxsta1(xx, yy, zz, vals[SW_VZ], sub->r1, sub->r2);
  case BX: return vals[SW_BX];
  case BY: return vals[SW_BY];
  case BZ: return vals[SW_BZ];
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_obstacle_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_obstacle, x)
static struct param ggcm_mhd_ic_obstacle_descr[] = {
  { "r1"           , VAR(r1)           , PARAM_FLOAT(1.)           },
  { "r2"           , VAR(r2)           , PARAM_FLOAT(2.)           },

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
// ggcm_mhd_ic subclass "obstacle"

struct ggcm_mhd_ic_ops ggcm_mhd_ic_obstacle_ops = {
  .name             = ggcm_mhd_ic_obstacle_name,
  .size             = sizeof(struct ggcm_mhd_ic_obstacle),
  .param_descr      = ggcm_mhd_ic_obstacle_descr,
  .primitive        = ggcm_mhd_ic_obstacle_primitive,
};
