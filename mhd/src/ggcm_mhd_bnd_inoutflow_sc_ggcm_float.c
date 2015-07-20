
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>
#include <mrc_fld_as_float.h>

#include <assert.h>

#define MT MT_SEMI_CONSERVATIVE_GGCM
#define SHIFT -1

#include "ggcm_mhd_bnd_inoutflow_common.c"

// ======================================================================
// ggcm_mhd_bnd subclass "inoutflow_sc_ggcm_float"

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_float = {
  .name             = "inoutflow_sc_ggcm_float",
  .fill_ghosts      = ggcm_mhd_bnd_sub_fill_ghosts,
};

