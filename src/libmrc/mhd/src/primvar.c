
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

#include <string.h>

// ----------------------------------------------------------------------
// primvar_c

void
primvar_c(struct ggcm_mhd *mhd, int m)
{
  const char *type = mrc_fld_type(mhd->fld);

  if (strcmp(type, "float") == 0) {
    return primvar_float(mhd, m);
  } else if (strcmp(type, "double") == 0) {
    return primvar_double(mhd, m);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// primvar1_c

void
primvar1_c(struct ggcm_mhd *mhd)
{
  return primvar_c(mhd, 0);
}

