
#include "psc_push_particles_ps.h"
#include <mrc_profile.h>

#define SIMD_BITS 0
#define PFX(x) sb0_ ## x

#include "ps_common.c"

#include "ps_1vb_yz.c"

