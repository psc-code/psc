
#include "ggcm_mhd_bndsw_private.h"

#include "ggcm_mhd_defs.h"

#include <assert.h>

// ======================================================================
// ggcm_mhd_bndsw class

// ----------------------------------------------------------------------
// ggcm_mhd_bndsw_new_step
//
// called at beginning of a new time step to reset cache pointer

void
ggcm_mhd_bndsw_new_step(struct ggcm_mhd_bndsw *bndsw)
{
  struct ggcm_mhd_bndsw_ops *ops = ggcm_mhd_bndsw_ops(bndsw);
  assert(ops);
  if (ops->new_step) {
    ops->new_step(bndsw);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bndsw_at

void
ggcm_mhd_bndsw_at(struct ggcm_mhd_bndsw *bndsw, float bntim, float xx[3],
		  float vals[])
{
  struct ggcm_mhd_bndsw_ops *ops = ggcm_mhd_bndsw_ops(bndsw);
  assert(ops && ops->at);
  ops->at(bndsw, bntim, xx, vals);
  vals[RR] *= bndsw->alphafak;
}

// ----------------------------------------------------------------------
// ggcm_mhd_bndsw_get_initial

void
ggcm_mhd_bndsw_get_initial(struct ggcm_mhd_bndsw *bndsw, float vals[])
{
  struct ggcm_mhd_bndsw_ops *ops = ggcm_mhd_bndsw_ops(bndsw);
  assert(ops && ops->get_initial);
  ops->get_initial(bndsw, vals);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bndsw_init

static void
ggcm_mhd_bndsw_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bndsw,
			      &ggcm_mhd_bndsw_none_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bndsw description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bndsw, x)
static struct param ggcm_mhd_bndsw_descr[] = {
  { "alphafak"        , VAR(alphafak)        , PARAM_FLOAT(1.)          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_bndsw class description

struct mrc_class_ggcm_mhd_bndsw mrc_class_ggcm_mhd_bndsw = {
  .name             = "ggcm_mhd_bndsw",
  .size             = sizeof(struct ggcm_mhd_bndsw),
  .param_descr      = ggcm_mhd_bndsw_descr,
  .init             = ggcm_mhd_bndsw_init,
};

