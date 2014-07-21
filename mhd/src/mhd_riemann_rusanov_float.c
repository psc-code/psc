
#include <mrc_fld_as_float.h>
#define F1(f, m, i) MRC_S2(f, m, i)

#include "mhd_riemann_rusanov_common.c"

// ----------------------------------------------------------------------
// mhd_riemann_rusanov_float_ops

struct mhd_riemann_ops mhd_riemann_rusanov_float_ops = {
  .name             = "rusanov_float",
  .run              = mhd_riemann_rusanov_run,
};

