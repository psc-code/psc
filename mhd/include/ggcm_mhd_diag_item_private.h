
#ifndef GGCM_MHD_DIAG_ITEM_PRIVATE_H
#define GGCM_MHD_DIAG_ITEM_PRIVATE_H

#include "ggcm_mhd_diag_item.h"

struct ggcm_mhd_diag_item {
  struct mrc_obj obj;
  struct ggcm_mhd_diag *diag;
};

struct ggcm_mhd_diag_item_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_diag_item);
  void (*run)(struct ggcm_mhd_diag_item *item,
	      struct mrc_io *io, struct mrc_fld *f,
	      int diag_type, float plane);
};

extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rr1;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_uu1;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_ee1;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rv1;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b1;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rr;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_v;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_pp;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_j;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_divb;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rank;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_e_ec;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_e_cc;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_gkeyll_e;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_gkeyll_i;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_gkeyll_em;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_ymask;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_zmask;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rmask;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_bnd_mask;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b0;
extern struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_psi;

#endif
