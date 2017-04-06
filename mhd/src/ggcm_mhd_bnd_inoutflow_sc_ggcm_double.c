
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>
#include <mrc_fld_as_double.h>

#include <assert.h>

#define MT MT_SCONS_FC_GGCM
#define SHIFT -1

#define ggcm_mhd_bnd_ops_inoutflow ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_double
#define ggcm_mhd_bnd_sub_name "inoutflow_sc_ggcm_double"

#include "ggcm_mhd_bnd_inoutflow_common.c"
