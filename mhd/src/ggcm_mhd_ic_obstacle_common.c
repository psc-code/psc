#include <stdlib.h>
#include <math.h>

#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bndsw.h"

#include "mrc_domain.h"

// there's no need to init B using vector potential, B being constant,
// but it can be useful for testing
#define USE_VECTOR_POTENTIAL 1

// ======================================================================
// ggcm_mhd_ic subclass "obstacle"

struct ggcm_mhd_ic_obstacle {
  float r1; // radius where velocity is ramped down all the way to 0
  float r2; // radius where velocity ramp down starts

  // solar wind values
  double bnvals[N_PRIMITIVE];
};

#define ggcm_mhd_ic_obstacle(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_obstacle)

// ----------------------------------------------------------------------
// get_solar_wind

static void
get_solar_wind(struct ggcm_mhd_ic *ic, mrc_fld_data_t vals[])
{
  struct ggcm_mhd_ic_obstacle *sub = ggcm_mhd_ic_obstacle(ic);
  struct ggcm_mhd *mhd = ic->mhd;

  vals[RR] = sub->bnvals[RR] / mhd->rrnorm;
  vals[VX] = sub->bnvals[VX] / mhd->vvnorm;
  vals[VY] = sub->bnvals[VY] / mhd->vvnorm;
  vals[VZ] = sub->bnvals[VZ] / mhd->vvnorm;
  vals[PP] = sub->bnvals[PP] / mhd->ppnorm;
  vals[BX] = sub->bnvals[BX] / mhd->bbnorm;
  vals[BY] = sub->bnvals[BY] / mhd->bbnorm;
  vals[BZ] = sub->bnvals[BZ] / mhd->bbnorm;
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

  mrc_fld_data_t vals[N_PRIMITIVE];
  get_solar_wind(ic, vals);

  double xx = crd[0], yy = crd[1], zz = crd[2];

  switch (m) {
  case RR: return vals[RR];
  case PP: return vals[PP];
  case VX: return vxsta1(xx, yy, zz, vals[VX], sub->r1, sub->r2);
  case VY: return vxsta1(xx, yy, zz, vals[VY], sub->r1, sub->r2);
  case VZ: return vxsta1(xx, yy, zz, vals[VZ], sub->r1, sub->r2);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_obstacle_primitive_bg

static double
ggcm_mhd_ic_obstacle_primitive_bg(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  mrc_fld_data_t vals[N_PRIMITIVE];
  get_solar_wind(ic, vals);

  switch (m) {
  case BX: return vals[BX];
  case BY: return vals[BY];
  case BZ: return vals[BZ];
  default: return 0.;
  }
}

#ifdef USE_VECTOR_POTENTIAL

// ----------------------------------------------------------------------
// ggcm_mhd_ic_obstacle_vector_potential_bg

static double
ggcm_mhd_ic_obstacle_vector_potential_bg(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  mrc_fld_data_t vals[N_PRIMITIVE];
  get_solar_wind(ic, vals);

  double xx = crd[0], yy = crd[1];

  switch (m) {
  case 1: return vals[BZ] * xx;
  case 2: return vals[BX] * yy - vals[BY] * xx;
  default: return 0.;
  }
}

#endif

// ----------------------------------------------------------------------
// ggcm_mhd_ic_obstacle_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_obstacle, x)
static struct param ggcm_mhd_ic_obstacle_descr[] = {
  { "r1"           , VAR(r1)           , PARAM_FLOAT(1.)           },
  { "r2"           , VAR(r2)           , PARAM_FLOAT(2.)           },

  { "rr"           , VAR(bnvals[RR])   , PARAM_DOUBLE(1.)          },
  { "pp"           , VAR(bnvals[PP])   , PARAM_DOUBLE(1.)          },
  { "vx"           , VAR(bnvals[VX])   , PARAM_DOUBLE(0.)          },
  { "vy"           , VAR(bnvals[VY])   , PARAM_DOUBLE(0.)          },
  { "vz"           , VAR(bnvals[VZ])   , PARAM_DOUBLE(0.)          },
  { "bx"           , VAR(bnvals[BX])   , PARAM_DOUBLE(0.)          },
  { "by"           , VAR(bnvals[BY])   , PARAM_DOUBLE(0.)          },
  { "bz"           , VAR(bnvals[BZ])   , PARAM_DOUBLE(0.)          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic subclass "obstacle"

struct ggcm_mhd_ic_ops ggcm_mhd_ic_obstacle_ops = {
  .name                = ggcm_mhd_ic_obstacle_name,
  .size                = sizeof(struct ggcm_mhd_ic_obstacle),
  .param_descr         = ggcm_mhd_ic_obstacle_descr,
  .primitive           = ggcm_mhd_ic_obstacle_primitive,
  .primitive_bg        = ggcm_mhd_ic_obstacle_primitive_bg,
#ifdef USE_VECTOR_POTENTIAL
  .vector_potential_bg = ggcm_mhd_ic_obstacle_vector_potential_bg,
#endif
};
