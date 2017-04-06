
#include "ggcm_mhd_step_private.h"

#include <ggcm_mhd_step.h>
#include "ggcm_mhd_step_gkeyll_lua.h"
#include <ggcm_mhd_gkeyll.h>

#include <ggcm_mhd_private.h>
#include <mrc_fld_as_double.h>

#include <string.h>

// ======================================================================
// ggcm_mhd_step subclass "gkeyll"

struct ggcm_mhd_step_gkeyll {
  const char *script;
  bool has_ymask;
  bool background;
  void *lua_state; // not using LuaState * b/c it is c++

  // intermediate states for dimension-splitting
  struct mrc_fld *qFlds[3];
};

#define ggcm_mhd_step_gkeyll(step) mrc_to_subobj(step, struct ggcm_mhd_step_gkeyll)

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_setup_flds

static void
ggcm_mhd_step_gkeyll_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_bool(mhd->fld, "aos", true);
  mrc_fld_set_param_bool(mhd->fld, "c_order", true);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_GKEYLL);

  ggcm_mhd_step_gkeyll_setup_flds_lua(mhd);

}

static void
ggcm_mhd_setup_gk_extra(struct ggcm_mhd *mhd) {
  ggcm_mhd_gkeyll_fluid_species_index_all(mhd, mhd->par.gk_idx);
  ggcm_mhd_gkeyll_fluid_species_q_m_all(mhd, mhd->par.gk_q_m);
  ggcm_mhd_gkeyll_fluid_species_mass_ratios_all(mhd, mhd->par.gk_mass_ratios);

  //re-normalize pressure ratios
  float pressure_total = 0.;
  for (int sp = 0; sp < mhd->par.gk_nr_fluids; sp++)
    pressure_total += mhd->par.gk_pressure_ratios.vals[sp];
  for (int sp = 0; sp < mhd->par.gk_nr_fluids; sp++)
    mhd->par.gk_pressure_ratios.vals[sp] /= pressure_total;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_setup

static void
ggcm_mhd_step_gkeyll_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_gkeyll *sub = ggcm_mhd_step_gkeyll(step);
  struct ggcm_mhd *mhd = step->mhd;

  ggcm_mhd_setup_gk_extra(mhd);

  if (sub->has_ymask) {
    mhd->ymask = ggcm_mhd_get_3d_fld(mhd, 1);
    mrc_fld_set(mhd->ymask, 1.);
  }
  if (sub->background) {
    mhd->b0 = ggcm_mhd_get_3d_fld(mhd, 3);
  }

  const int *dims = mrc_fld_spatial_dims(mhd->fld);
  sub->qFlds[0] = ggcm_mhd_get_3d_fld(mhd, mhd->fld->_nr_comps);
  mrc_fld_dict_add_int(sub->qFlds[0], "mhd_type", MT_GKEYLL);
  if (dims[1] > 1) {
    sub->qFlds[1] = ggcm_mhd_get_3d_fld(mhd, mhd->fld->_nr_comps);
    mrc_fld_dict_add_int(sub->qFlds[1], "mhd_type", MT_GKEYLL);
  }
  if (dims[2] > 1) {
    sub->qFlds[2] = ggcm_mhd_get_3d_fld(mhd, mhd->fld->_nr_comps);
    mrc_fld_dict_add_int(sub->qFlds[2], "mhd_type", MT_GKEYLL);
  }

  assert(strcmp(mrc_fld_type(mhd->fld), FLD_TYPE) == 0);
  assert(mhd->fld->_aos);
  assert(mhd->fld->_c_order);
  ggcm_mhd_step_gkeyll_lua_setup(&(sub->lua_state), sub->script, mhd, sub->qFlds);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_destroy

static void
ggcm_mhd_step_gkeyll_destroy(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_gkeyll *sub = ggcm_mhd_step_gkeyll(step);
  struct ggcm_mhd *mhd = step->mhd;

  ggcm_mhd_step_gkeyll_lua_destroy(sub->lua_state, mhd);

  if (sub->has_ymask) {
    ggcm_mhd_put_3d_fld(mhd, mhd->ymask);
  }
  if (sub->background) {
    ggcm_mhd_put_3d_fld(mhd, mhd->b0);
  }

  for (int d = 0; d < 3; d++)
     ggcm_mhd_put_3d_fld(mhd, sub->qFlds[d]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll_run

static void
ggcm_mhd_step_gkeyll_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_gkeyll *sub = ggcm_mhd_step_gkeyll(step);
  struct ggcm_mhd *mhd = step->mhd;

  assert(strcmp(mrc_fld_type(x), "double") == 0);
  assert(x->_aos && x->_c_order);
  ggcm_mhd_step_gkeyll_lua_run(sub->lua_state, mhd, x);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_gkeyll subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_gkeyll, x)
static struct param ggcm_mhd_step_gkeyll_descr[] = {
  { "script"       , VAR(script),        PARAM_STRING("gkeyll_step.lua") },
  { "has_ymask"    , VAR(has_ymask),     PARAM_BOOL(false)               },
  { "background"   , VAR(background),    PARAM_BOOL(false)               },
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

