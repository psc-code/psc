
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

#include <mrc_fld_as_double.h>

#define MT MT_BGRID_FC
#define ggcm_mhd_calc_currcc_bgrid ggcm_mhd_calc_currcc_bgrid_fc

#include "calc_currcc_bgrid_common.c"
