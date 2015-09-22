
#ifndef GGCM_MHD_BND_PRIVATE_H
#define GGCM_MHD_BND_PRIVATE_H

#include "ggcm_mhd_bnd.h"

struct ggcm_mhd_bnd {
  struct mrc_obj obj;
  struct ggcm_mhd *mhd;
};

struct ggcm_mhd_bnd_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_bnd);
  void (*fill_ghosts)(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
		      int m, float bntim);
};

extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_none;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_fc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_double_aos;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_fc_double_aos;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_gkeyll;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_sc_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_fc_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_sc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_fc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_sc_double_aos;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_fc_double_aos;

#define ggcm_mhd_bnd_ops(bnd) ((struct ggcm_mhd_bnd_ops *)((bnd)->obj.ops))

#endif
