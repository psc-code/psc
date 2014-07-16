
#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include "mhd_reconstruct_pcm_common.c"

// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_ops

struct mhd_reconstruct_ops mhd_reconstruct_pcm_double_ops = {
  .name             = "pcm_double",
  .run              = mhd_reconstruct_pcm_run,
};

