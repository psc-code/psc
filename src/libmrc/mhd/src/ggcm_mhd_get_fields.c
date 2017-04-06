#include "ggcm_mhd_step_cweno_private.h" 
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_io.h>

// ----------------------------------------------------------------------
// ggcm_mhd_get_fields

struct mrc_fld *
ggcm_mhd_get_fields(struct ggcm_mhd *mhd, const char *name, int nr_comps)
{   
  struct mrc_fld *f3 = mrc_domain_fld_create(mhd->domain, SW_2, NULL);
  mrc_fld_set_name(f3, name);
  mrc_fld_set_param_int(f3, "nr_comps", nr_comps);
  mrc_fld_setup(f3);
  return f3;
}
