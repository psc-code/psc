
#include <mrc_fld_as_double.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include "ggcm_mhd_step_c3_common.c"

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c3_double"

struct ggcm_mhd_step_ops ggcm_mhd_step_c3_double_ops = {
  .name        = "c3_double",
  .size        = sizeof(struct ggcm_mhd_step_c3),
  .param_descr = ggcm_mhd_step_c_descr,
  .mhd_type    = MT_SEMI_CONSERVATIVE,
  .fld_type    = FLD_TYPE,
  .nr_ghosts   = 2,
  .create      = ggcm_mhd_step_c_create,
  .setup       = ggcm_mhd_step_c_setup,
  .destroy     = ggcm_mhd_step_c_destroy,
  .run         = ggcm_mhd_step_c_run,
  .get_e_ec   = ggcm_mhd_step_c3_get_e_ec,
};
