
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_crds.h"

#include <mrc_profile.h>
#include <mrc_io.h>

#include <math.h>
#include <string.h>
#include <assert.h>

#include <mrc_fld_as_float.h>

#include "mhd_sc.c"

void
newstep(struct ggcm_mhd *mhd, float *dtn)
{
  ggcm_mhd_fill_ghosts(mhd, mhd->fld, _RR1, mhd->time);

  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);
  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM) {
    newstep_c2(mhd, dtn);
  } else if (mhd_type == MT_SEMI_CONSERVATIVE) {
    zmaskn(mhd, mhd->fld);
    if (strcmp(mrc_fld_type(mhd->fld), "double") == 0) {
      assert(0);  // USE step_c3_double instead
    }
    *dtn = newstep_sc(mhd, mhd->fld);
  } else {
    assert(0);
  }
}
