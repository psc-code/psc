
#include <mrc_fld_as_float.h>

#include "primvar_common.c"

// ----------------------------------------------------------------------
// primvar1_c

void
primvar1_c(struct ggcm_mhd *mhd)
{
  return primvar_c(mhd, _RR1);
}

