
#include <mrc_domain.h>
#include <mrc_profile.h>
#include <string.h>

static inline bool
at_lo_boundary(int p, int d)
{
  return ppsc->patch[p].off[d] == 0;
}

static inline bool
at_hi_boundary(int p, int d)
{
  return ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d];
}

#define DDCP_TYPE DDCP_TYPE_COMMON_OMP

#include "ddc_particles_inc.c"

