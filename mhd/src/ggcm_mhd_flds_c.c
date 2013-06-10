
#include "ggcm_mhd_flds_private.h"

#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_io.h>
#include <assert.h>
#include <string.h>

// ======================================================================
// ggcm_mhd_flds subclass "c"

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_create

static void
ggcm_mhd_flds_c_create(struct ggcm_mhd_flds *flds)
{
  mrc_fld_set_type(flds->fld, "float");
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds subclass "c"

struct ggcm_mhd_flds_ops ggcm_mhd_flds_ops_c = {
  .name             = "c",
  .create           = ggcm_mhd_flds_c_create,
};
