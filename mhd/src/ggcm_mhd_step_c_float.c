
#include <mrc_fld_as_float.h>

#include "ggcm_mhd_step_c_common.c"

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c_float"

struct ggcm_mhd_step_ops ggcm_mhd_step_c_float_ops = {
  .name        = "c_float",
  .mhd_type    = MT_SEMI_CONSERVATIVE_GGCM,
  .fld_type    = FLD_TYPE,
  .nr_ghosts   = 2,
  .newstep     = ggcm_mhd_step_c_newstep,
  .pred        = ggcm_mhd_step_c_pred,
  .corr        = ggcm_mhd_step_c_corr,
  .run         = ggcm_mhd_step_run_predcorr,
  .get_e_ec   = ggcm_mhd_step_c_get_e_ec,
};
