
#include <mrc_fld_as_float.h>
#define _F1(f, m, i) MRC_S2(f, m, i)

#define ggcm_mhd_step_c3_ops ggcm_mhd_step_c3_float_ops
#define ggcm_mhd_step_c3_name "c3_float"

#include "ggcm_mhd_step_c3_common.c"
