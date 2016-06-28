
#include <mrc_fld_as_float.h>

#define ggcm_mhd_step_c_ops ggcm_mhd_step_c_float_ops
#define ggcm_mhd_step_c_name "c_float"

#define OPT_STAGGER OPT_STAGGER_GGCM

#include "ggcm_mhd_step_c_common.c"

