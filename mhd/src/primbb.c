
#include "ggcm_mhd_private.h"

void
primbb_c(struct ggcm_mhd *mhd, int m_curr)
{
  return primbb_float(mhd, m_curr);
}

void
primbb_c2_c(struct ggcm_mhd *mhd, int m_curr)
{
  return primbb_c2_float(mhd, m_curr);
}
