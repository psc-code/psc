
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>
#include <mrc_fld_as_float.h>

#include <assert.h>

#define MT MT_SCONS_FC_GGCM
#define SHIFT -1

#define ggcm_mhd_bnd_ops_inoutflow ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_float
#define ggcm_mhd_bnd_sub_name "inoutflow_sc_ggcm_float"

#include "ggcm_mhd_bnd_inoutflow_common.c"
