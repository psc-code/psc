
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>
#include <mrc_fld_as_double.h>

#include <assert.h>

#define MT MT_SEMI_CONSERVATIVE
#define SHIFT 0

#include "ggcm_mhd_bnd_inoutflow_common.c"

// ======================================================================
// ggcm_mhd_bnd subclass "inoutflow_sc_double"

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_double = {
  .name             = "inoutflow_sc_double",
  .fill_ghosts      = ggcm_mhd_bnd_sub_fill_ghosts,
};

