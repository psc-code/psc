
#include <mrc_fld_as_float.h>

#define copy_sc_ggcm_to_sc copy_sc_ggcm_to_sc_float
#define copy_sc_ggcm_to_fc copy_sc_ggcm_to_fc_float
#define copy_sc_to_sc_ggcm copy_sc_to_sc_ggcm_float
#define copy_fc_to_sc_ggcm copy_fc_to_sc_ggcm_float
#define copy_fc_to_sc copy_fc_to_sc_float
#define copy_sc_to_fc copy_sc_to_fc_float

#include "ggcm_mhd_convert_common.c"
