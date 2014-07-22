
#include <mrc_fld_as_float.h>
#define F1(f, m, i) MRC_S2(f, m, i)

#include "mhd_reconstruct_plm_common.c"

// ----------------------------------------------------------------------
// mhd_reconstruct_plm_ops

struct mhd_reconstruct_ops mhd_reconstruct_plm_float_ops = {
  .name             = "plm_float",
  .run              = mhd_reconstruct_plm_run,
};

