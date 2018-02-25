
#include "psc_generic_c.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct Config2ndZ
{
  using mparticles_t = PscMparticlesDouble;
  using mfields_t = PscMfieldsC;
};

#define CONFIG Config2ndZ

#define DIM DIM_Z
#define ORDER ORDER_2ND
#define PRTS PRTS_STAGGERED
#include "push_part_common.c"

