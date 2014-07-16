
#include "ggcm_mhd_private.h"

#include <string.h>

void
zmaskn_c(struct ggcm_mhd *mhd)
{
  const char *type = mrc_fld_type(mhd->fld);

  if (strcmp(type, "float") == 0) {
    return zmaskn_float(mhd);
  } else if (strcmp(type, "double") == 0) {
    return zmaskn_double(mhd);
  } else {
    assert(0);
  }
}

