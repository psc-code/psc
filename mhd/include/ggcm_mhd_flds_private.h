
#ifndef GGCM_MHD_FLDS_PRIVATE_H
#define GGCM_MHD_FLDS_PRIVATE_H

#include "ggcm_mhd_flds.h"

#include "ggcm_mhd_defs.h"

struct ggcm_mhd_flds {
  struct mrc_obj obj;
  struct mrc_fld *fld;
};

struct ggcm_mhd_flds_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_flds);
  void (*copy)(struct ggcm_mhd_flds *to, struct ggcm_mhd_flds *from);
};

typedef void (*ggcm_mhd_flds_copy_to_func_t)(struct ggcm_mhd_flds *,
					     struct ggcm_mhd_flds *);
typedef void (*ggcm_mhd_flds_copy_from_func_t)(struct ggcm_mhd_flds *,
					       struct ggcm_mhd_flds *);

#define ggcm_mhd_flds_ops(flds) ((struct ggcm_mhd_flds_ops *)(flds)->obj.ops)

extern struct ggcm_mhd_flds_ops ggcm_mhd_flds_ops_c;

// ======================================================================
// this is the subclass struct for both aos and naos

struct ggcm_mhd_flds_cgen {
  float *flds;
  float *_flds1, *_flds2, *_fldsp, *_dummy_dma;
};

#define ggcm_mhd_flds_cgen(flds) mrc_to_subobj(flds, struct ggcm_mhd_flds_cgen)

#endif


