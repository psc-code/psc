
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>
#include <mrc_fld_as_double.h>

#include <assert.h>

#define MT MT_GKEYLL
#define SHIFT 0

#define ggcm_mhd_bnd_ops_inoutflow ggcm_mhd_bnd_ops_inoutflow_gkeyll
#define ggcm_mhd_bnd_sub_name "inoutflow_gkeyll"

#include "ggcm_mhd_bnd_inoutflow_common.c"
