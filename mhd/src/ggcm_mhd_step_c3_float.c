#define BOUNDS_CHECK
#include <mrc_fld_as_float.h>
#define F1(f, m, i) MRC_S2(f, m, i)

#include "ggcm_mhd_step_c3_common.c"

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c3_float"

struct ggcm_mhd_step_ops ggcm_mhd_step_c3_float_ops = {
  .name        = "c3_float",
  .size        = sizeof(struct ggcm_mhd_step_c3),
  .param_descr = ggcm_mhd_step_c_descr,
  .mhd_type    = MT_SEMI_CONSERVATIVE,
  .fld_type    = FLD_TYPE,
  .nr_ghosts   = 2,
  .setup       = ggcm_mhd_step_c_setup,
  .run         = ggcm_mhd_step_c_run,
};
