
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

// ----------------------------------------------------------------------
// primvar_c

void
primvar_c(struct ggcm_mhd *mhd, int m)
{
  return primvar_float(mhd, m);
}

// ----------------------------------------------------------------------
// primvar1_c

void
primvar1_c(struct ggcm_mhd *mhd)
{
  return primvar_c(mhd, _RR1);
}

