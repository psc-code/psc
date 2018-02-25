
#include "psc_generic_c.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct Config2ndXZ
{
  using mparticles_t = PscMparticlesDouble;
  using mfields_t = PscMfieldsC;
};

#define CONFIG Config2ndXZ

#define DIM DIM_XZ
#define ORDER ORDER_2ND
#define PRTS PRTS_STAGGERED
#include "push_part_common.c"

