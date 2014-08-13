
#ifndef MHD_RIEMANN_PRIVATE_H
#define MHD_RIEMANN_PRIVATE_H

#include "mhd_riemann.h"

struct mhd_riemann {
  struct mrc_obj obj;
  struct ggcm_mhd *mhd;
};

struct mhd_riemann_ops {
  MRC_SUBCLASS_OPS(struct mhd_riemann);
  void (*run)(struct mhd_riemann *riemann,
	      struct mrc_fld *F1d,
	      struct mrc_fld *Ul, struct mrc_fld *Ur,
	      struct mrc_fld *Wl, struct mrc_fld *Wr,
	      int ldim, int l, int r, int dir);
};

#define mhd_riemann_ops(ai) ((struct mhd_riemann_ops *)(ai)->obj.ops)

extern struct mhd_riemann_ops mhd_riemann_rusanov_double_ops;
extern struct mhd_riemann_ops mhd_riemann_rusanov_float_ops;
extern struct mhd_riemann_ops mhd_riemann_hll_ops;
extern struct mhd_riemann_ops mhd_riemann_hlld_ops;
extern struct mhd_riemann_ops mhd_riemann_hydro_rusanov_ops;
extern struct mhd_riemann_ops mhd_riemann_hydro_hll_ops;
extern struct mhd_riemann_ops mhd_riemann_hydro_hllc_ops;

#endif
