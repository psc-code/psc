
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd.h"

#include <mrc_io.h>

#include <assert.h>

// ======================================================================
// ggcm_mhd_bnd class

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_fill_ghosts

void
ggcm_mhd_bnd_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
			 float bntim)
{
  struct ggcm_mhd_bnd_ops *ops = ggcm_mhd_bnd_ops(bnd);
  assert(ops && ops->fill_ghosts);
  ops->fill_ghosts(bnd, fld, bntim);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_fill_ghosts_E

void
ggcm_mhd_bnd_fill_ghosts_E(struct ggcm_mhd_bnd *bnd, struct mrc_fld *E)
{
  struct ggcm_mhd_bnd_ops *ops = ggcm_mhd_bnd_ops(bnd);
  assert(ops);
  if (ops->fill_ghosts_E) {
    ops->fill_ghosts_E(bnd, E);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_fill_ghosts_reconstr

void
ggcm_mhd_bnd_fill_ghosts_reconstr(struct ggcm_mhd_bnd *bnd, struct mrc_fld *U_l[],
				  struct mrc_fld *U_r[], int p)
{
  struct ggcm_mhd_bnd_ops *ops = ggcm_mhd_bnd_ops(bnd);
  assert(ops);
  if (ops->fill_ghosts_reconstr) {
    ops->fill_ghosts_reconstr(bnd, U_l, U_r, p);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_init

static void
ggcm_mhd_bnd_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_none);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_float);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_inoutflow_sc_float);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_inoutflow_sc_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_inoutflow_fc_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_inoutflow_fc_cc_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_inoutflow_gkeyll);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_sphere_sc_float);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_sphere_fc_float);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_sphere_sc_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_sphere_sc_ggcm_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_sphere_fc_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_sphere_fc_cc_double);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_ops_sphere_gkeyll);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bnd, x)
static struct param ggcm_mhd_bnd_descr[] = {
  { "mhd"             , VAR(mhd)             , PARAM_OBJ(ggcm_mhd)      },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_bnd class description

struct mrc_class_ggcm_mhd_bnd mrc_class_ggcm_mhd_bnd = {
  .name             = "ggcm_mhd_bnd",
  .size             = sizeof(struct ggcm_mhd_bnd),
  .param_descr      = ggcm_mhd_bnd_descr,
  .init             = ggcm_mhd_bnd_init,
};

