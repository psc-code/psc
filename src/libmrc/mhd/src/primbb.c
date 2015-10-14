
#include "ggcm_mhd_private.h"

#include <string.h>

void
primbb_c(struct ggcm_mhd *mhd, int m)
{
  const char *type = mrc_fld_type(mhd->fld);

  if (strcmp(type, "float") == 0) {
    return primbb_float(mhd, m);
  } else if (strcmp(type, "double") == 0) {
    return primbb_double(mhd, m);
  } else {
    assert(0);
  }
}
