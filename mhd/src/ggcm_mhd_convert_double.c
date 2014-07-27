
#include <mrc_fld_as_double.h>

#define copy_sc_ggcm_to_sc copy_sc_ggcm_to_sc_double
#define copy_sc_ggcm_to_fc copy_sc_ggcm_to_fc_double
#define copy_sc_to_sc_ggcm copy_sc_to_sc_ggcm_double
#define copy_fc_to_sc_ggcm copy_fc_to_sc_ggcm_double
#define copy_fc_to_sc copy_fc_to_sc_double
#define copy_sc_to_fc copy_sc_to_fc_double

#include "ggcm_mhd_convert_common.c"
