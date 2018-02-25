
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct Config1stYZ
{
  using mparticles_t = PscMparticlesDouble;
  using mfields_t = PscMfieldsC;
};

#define CONFIG Config1stYZ

#define DIM DIM_YZ
#define ORDER ORDER_1ST
#define PRTS PRTS_STAGGERED
#include "push_part_common.c"

