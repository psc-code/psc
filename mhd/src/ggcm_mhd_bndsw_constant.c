
#include "ggcm_mhd_bndsw_private.h"

struct ggcm_mhd_bndsw_constant {
  double bnvals[SW_NR];
};

static void ggcm_mhd_bndsw_constant_at(struct ggcm_mhd_bndsw *bndsw, 
  float bntim, float xx[3], float vals[SW_NR]) {
  struct ggcm_mhd_bndsw_constant *sub = mrc_to_subobj(bndsw, struct ggcm_mhd_bndsw_constant);
  vals[SW_RR] = sub->bnvals[SW_RR];
  vals[SW_VX] = sub->bnvals[SW_VX];
  vals[SW_VY] = sub->bnvals[SW_VY];
  vals[SW_VZ] = sub->bnvals[SW_VZ];
  vals[SW_PP] = sub->bnvals[SW_PP];
  vals[SW_BX] = sub->bnvals[SW_BX];
  vals[SW_BY] = sub->bnvals[SW_BY];
  vals[SW_BZ] = sub->bnvals[SW_BZ];
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd inflow description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bndsw_constant, x)
static struct param ggcm_mhd_bndsw_constant_descr[] = {
  { "rr"           , VAR(bnvals[SW_RR])           , PARAM_DOUBLE(1.) },
  { "pp"           , VAR(bnvals[SW_PP])           , PARAM_DOUBLE(1.) },
  { "vx"           , VAR(bnvals[SW_VX])           , PARAM_DOUBLE(0.) },
  { "vy"           , VAR(bnvals[SW_VY])           , PARAM_DOUBLE(0.) },
  { "vz"           , VAR(bnvals[SW_VZ])           , PARAM_DOUBLE(0.) },
  { "bx"           , VAR(bnvals[SW_BX])           , PARAM_DOUBLE(0.) },
  { "by"           , VAR(bnvals[SW_BY])           , PARAM_DOUBLE(0.) },
  { "bz"           , VAR(bnvals[SW_BZ])           , PARAM_DOUBLE(0.) },

  {},
};
#undef VAR


// ----------------------------------------------------------------------
// ggcm_mhd_bndsw subclass "constant"

struct ggcm_mhd_bndsw_ops ggcm_mhd_bndsw_constant_ops = {
  .name             = "constant",
  .size             = sizeof(struct ggcm_mhd_bndsw_constant),
  .param_descr      = ggcm_mhd_bndsw_constant_descr,
  .at               = ggcm_mhd_bndsw_constant_at,
};

