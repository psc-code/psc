
#include <mrc_fld_as_float.h>

#include "ggcm_mhd_step_c3_common.c"

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c3_float"

struct ggcm_mhd_step_ops ggcm_mhd_step_c3_float_ops = {
  .name        = "c3_float",
  .mhd_type    = MT_SEMI_CONSERVATIVE,
  .fld_type    = FLD_TYPE,
  .nr_ghosts   = 2,
  .pred        = ggcm_mhd_step_c_pred,
  .corr        = ggcm_mhd_step_c_corr,
  .run         = ggcm_mhd_step_run_predcorr,
};
