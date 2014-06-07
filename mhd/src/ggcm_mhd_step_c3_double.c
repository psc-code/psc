
#include <mrc_fld_as_double.h>

#include "ggcm_mhd_step_c3_common.c"

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c3_double"

struct ggcm_mhd_step_ops ggcm_mhd_step_c3_double_ops = {
  .name        = "c3_double",
  .mhd_type    = MT_SEMI_CONSERVATIVE,
  .fld_type    = FLD_TYPE,
  .nr_ghosts   = 2,
  .run         = ggcm_mhd_step_c_run,
};
