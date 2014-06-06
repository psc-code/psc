
#include <mrc_fld_as_double.h>

#include "ggcm_mhd_step_c2_common.c"

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c2_double"

struct ggcm_mhd_step_ops ggcm_mhd_step_c2_double_ops = {
  .name        = "c2_double",
  .mhd_type    = MT_SEMI_CONSERVATIVE,
  .fld_type    = FLD_TYPE,
  .nr_ghosts   = 2,
  .pred        = ggcm_mhd_step_c_pred,
  .corr        = ggcm_mhd_step_c_corr,
  .run         = ggcm_mhd_step_run_predcorr,
};
