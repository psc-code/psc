
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_ic_gkeyll_lua.h"

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>
#include <mrc_domain.h>
#include <mrc_fld_as_double_aos.h>

#include <math.h>

// ======================================================================
// ggcm_mhd_ic subclass "gkeyll"

struct ggcm_mhd_ic_gkeyll {
  const char *script;
  const char *script_common;
};

#define ggcm_mhd_ic_gkeyll(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_gkeyll)

// ----------------------------------------------------------------------
// ggcm_mhd_ic_gkeyll_run

static void
ggcm_mhd_ic_gkeyll_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_gkeyll *sub = ggcm_mhd_ic_gkeyll(ic);
  struct ggcm_mhd *mhd = ic->mhd;

  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  ggcm_mhd_ic_gkeyll_lua_run(sub->script, sub->script_common, mhd, fld);
  mrc_fld_put_as(fld, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_gkeyll subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_gkeyll, x)
static struct param ggcm_mhd_ic_gkeyll_descr[] = {
  { "script"           , VAR(script)           , PARAM_STRING("init.lua")   },
  { "script_common"    , VAR(script_common)    , PARAM_STRING("common.lua") },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_gkeyll_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_gkeyll_ops = {
  .name        = "gkeyll",
  .size        = sizeof(struct ggcm_mhd_ic_gkeyll),
  .param_descr = ggcm_mhd_ic_gkeyll_descr,
  .run         = ggcm_mhd_ic_gkeyll_run,
};

