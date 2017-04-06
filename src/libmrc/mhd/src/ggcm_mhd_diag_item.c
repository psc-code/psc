
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_diag.h"

#include <stdio.h>
#include <assert.h>

#define ggcm_mhd_diag_item_ops(item) ((struct ggcm_mhd_diag_item_ops *) (item)->obj.ops)

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_run

void
ggcm_mhd_diag_item_run(struct ggcm_mhd_diag_item *item, struct mrc_io *io,
		       struct mrc_fld *f, int diag_type, float plane)
{
  struct ggcm_mhd_diag_item_ops *ops = ggcm_mhd_diag_item_ops(item);
  assert(ops && ops->run);
  ops->run(item, io, f, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_init

static void
ggcm_mhd_diag_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rr1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_uu1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_ee1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rv1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_b1);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rr);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_v);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_pp);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_b);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_divb);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rank);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_e_ec);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_e_cc);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_j);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_e);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_i);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_em);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_ymask);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_zmask);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_bnd_mask);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_rmask);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_b0);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_psi);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_diag_item, x)
static struct param ggcm_mhd_diag_item_descr[] = {
  { "diag"            , VAR(diag)            , PARAM_OBJ(ggcm_mhd_diag) },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item class

struct mrc_class_ggcm_mhd_diag_item mrc_class_ggcm_mhd_diag_item = {
  .name             = "ggcm_mhd_diag_item",
  .size             = sizeof(struct ggcm_mhd_diag_item),
  .param_descr      = ggcm_mhd_diag_item_descr,
  .init             = ggcm_mhd_diag_init,
};

