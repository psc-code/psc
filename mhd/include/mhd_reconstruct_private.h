
#ifndef MHD_RECONSTRUCT_PRIVATE_H
#define MHD_RECONSTRUCT_PRIVATE_H

#include "mhd_reconstruct.h"

struct mhd_reconstruct {
  struct mrc_obj obj;
  struct ggcm_mhd *mhd;
};

struct mhd_reconstruct_ops {
  MRC_SUBCLASS_OPS(struct mhd_reconstruct);
  void (*run)(struct mhd_reconstruct *mr,
	      struct mrc_fld *Ul, struct mrc_fld *Ur,
	      struct mrc_fld *Wl, struct mrc_fld *Wr,
	      struct mrc_fld *W1d, struct mrc_fld *Bxi,
	      int ldim, int l, int r, int dir);
};

#define mhd_reconstruct_ops(ai) ((struct mhd_reconstruct_ops *)(ai)->obj.ops)

extern struct mhd_reconstruct_ops mhd_reconstruct_pcm_double_ops;
extern struct mhd_reconstruct_ops mhd_reconstruct_pcm_float_ops;
extern struct mhd_reconstruct_ops mhd_reconstruct_plm_double_ops;
extern struct mhd_reconstruct_ops mhd_reconstruct_plm_float_ops;

#endif
