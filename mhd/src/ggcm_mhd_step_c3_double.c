
#include <mrc_fld_as_double.h>
#define _F1(f, m, i) MRC_D2(f, m, i)

#define ggcm_mhd_step_c3_ops ggcm_mhd_step_c3_double_ops
#define ggcm_mhd_step_c3_name "c3_double"

#include "ggcm_mhd_step_c3_common.c"

