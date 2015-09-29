
#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include "mhd_riemann_rusanov_common.c"

// ----------------------------------------------------------------------
// mhd_riemann_rusanov_double_ops

struct mhd_riemann_ops mhd_riemann_rusanov_double_ops = {
  .name             = "rusanov_double",
  .run              = mhd_riemann_rusanov_run,
};

