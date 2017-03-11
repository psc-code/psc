
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

#include <mrc_fld_as_double.h>

#define MT MT_BGRID_FC_GGCM
#define BGRID_SFX(x) x ## _bgrid_fc_ggcm

#include "calc_currcc_bgrid_common.c"
