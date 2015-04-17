
#include "ggcm_mhd_step_private.h"

#include "ggcm_mhd_step_gkeyll_lua.h"

#include <ggcm_mhd_private.h>
#include <mrc_fld_as_double_aos.h>

#include <string.h>

// ======================================================================
// ggcm_mhd_step subclass "gkeyll"

struct ggcm_mhd_step_gkeyll {
  const char *script;
};

#define ggcm_mhd_step_gkeyll(step) mrc_to_subobj(step, struct ggcm_mhd_step_gkeyll)

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_setup_flds

static void
ggcm_mhd_step_gkeyll_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  // FIXME, not nice to hardcode 2 here, should be gotten from the script
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_FULLY_CONSERVATIVE);
  // FIXME, not nice to hardcode 18 here, should be gotten from the script
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 18);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_setup

static void
ggcm_mhd_step_gkeyll_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_gkeyll *sub = ggcm_mhd_step_gkeyll(step);
  struct ggcm_mhd *mhd = step->mhd;

  assert(strcmp(mrc_fld_type(mhd->fld), FLD_TYPE) == 0);
  ggcm_mhd_step_gkeyll_lua_setup(sub->script, mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_destroy

static void
ggcm_mhd_step_gkeyll_destroy(struct ggcm_mhd_step *step)
{
  MHERE; // TBD
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_run

static void
ggcm_mhd_step_gkeyll_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;

  assert(strcmp(mrc_fld_type(x), "double_aos") == 0);
  ggcm_mhd_step_gkeyll_lua_run(mhd, x);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_gkeyll, x)
static struct param ggcm_mhd_step_gkeyll_descr[] = {
  { "script"           , VAR(script)           , PARAM_STRING("step.lua")   },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_ops

struct ggcm_mhd_step_ops ggcm_mhd_step_gkeyll_ops = {
  .name             = "gkeyll",
  .size             = sizeof(struct ggcm_mhd_step_gkeyll),
  .param_descr      = ggcm_mhd_step_gkeyll_descr,
  .setup            = ggcm_mhd_step_gkeyll_setup,
  .destroy          = ggcm_mhd_step_gkeyll_destroy,
  .run              = ggcm_mhd_step_gkeyll_run,
  .setup_flds       = ggcm_mhd_step_gkeyll_setup_flds,
};

