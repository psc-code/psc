
#include <mrc_fld_as_float.h>
#define F1(f, m, i) MRC_S2(f, m, i)

#include "mhd_reconstruct_pcm_common.c"

// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_ops

struct mhd_reconstruct_ops mhd_reconstruct_pcm_float_ops = {
  .name             = "pcm_float",
  .run              = mhd_reconstruct_pcm_run,
};

