
#include "ggcm_mhd_bndsw_constant_5m.h"

struct ggcm_mhd_bndsw_constant_5m {
  double bnvals[SW5_NR];
};

static void ggcm_mhd_bndsw_constant_5m_at(struct ggcm_mhd_bndsw *bndsw, 
  float bntim, float xx[3], float vals[SW5_NR]) {
  struct ggcm_mhd_bndsw_constant_5m *sub = mrc_to_subobj(bndsw, struct ggcm_mhd_bndsw_constant_5m);
  vals[SW5_RRE] = sub->bnvals[SW5_RRE];
  vals[SW5_PPE] = sub->bnvals[SW5_PPE];
  vals[SW5_VXE] = sub->bnvals[SW5_VXE];
  vals[SW5_VYE] = sub->bnvals[SW5_VYE];
  vals[SW5_VZE] = sub->bnvals[SW5_VZE];
  vals[SW5_RRI] = sub->bnvals[SW5_RRI];
  vals[SW5_PPI] = sub->bnvals[SW5_PPI];
  vals[SW5_VXI] = sub->bnvals[SW5_VXI];
  vals[SW5_VYI] = sub->bnvals[SW5_VYI];
  vals[SW5_VZI] = sub->bnvals[SW5_VZI];
  vals[SW5_EX] = sub->bnvals[SW5_EX];
  vals[SW5_EY] = sub->bnvals[SW5_EY];
  vals[SW5_EZ] = sub->bnvals[SW5_EZ];
  vals[SW5_BX] = sub->bnvals[SW5_BX];
  vals[SW5_BY] = sub->bnvals[SW5_BY];
  vals[SW5_BZ] = sub->bnvals[SW5_BZ];
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd inflow description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bndsw_constant_5m, x)
static struct param ggcm_mhd_bndsw_constant_5m_descr[] = {
  { "rre"           , VAR(bnvals[SW5_RRE])           , PARAM_DOUBLE(1.) },
  { "ppe"           , VAR(bnvals[SW5_PPE])           , PARAM_DOUBLE(1.) },
  { "vxe"           , VAR(bnvals[SW5_VXE])           , PARAM_DOUBLE(0.) },
  { "vye"           , VAR(bnvals[SW5_VYE])           , PARAM_DOUBLE(0.) },
  { "vze"           , VAR(bnvals[SW5_VZE])           , PARAM_DOUBLE(0.) },
  { "rri"           , VAR(bnvals[SW5_RRI])           , PARAM_DOUBLE(1.) },
  { "ppe"           , VAR(bnvals[SW5_PPI])           , PARAM_DOUBLE(1.) },
  { "vxi"           , VAR(bnvals[SW5_VXI])           , PARAM_DOUBLE(0.) },
  { "vyi"           , VAR(bnvals[SW5_VYI])           , PARAM_DOUBLE(0.) },
  { "vzi"           , VAR(bnvals[SW5_VZI])           , PARAM_DOUBLE(0.) },
  { "ex"            , VAR(bnvals[SW5_EX])            , PARAM_DOUBLE(0.) },
  { "ey"            , VAR(bnvals[SW5_EY])            , PARAM_DOUBLE(0.) },
  { "ez"            , VAR(bnvals[SW5_EZ])            , PARAM_DOUBLE(0.) },
  { "bx"            , VAR(bnvals[SW5_BX])            , PARAM_DOUBLE(0.) },
  { "by"            , VAR(bnvals[SW5_BY])            , PARAM_DOUBLE(0.) },
  { "bz"            , VAR(bnvals[SW5_BZ])            , PARAM_DOUBLE(0.) },

  {},
};
#undef VAR


// ----------------------------------------------------------------------
// ggcm_mhd_bndsw subclass "constant_5m"

struct ggcm_mhd_bndsw_ops ggcm_mhd_bndsw_constant_5m_ops = {
  .name             = "constant_5m",
  .size             = sizeof(struct ggcm_mhd_bndsw_constant_5m),
  .param_descr      = ggcm_mhd_bndsw_constant_5m_descr,
  .at               = ggcm_mhd_bndsw_constant_5m_at,
};

