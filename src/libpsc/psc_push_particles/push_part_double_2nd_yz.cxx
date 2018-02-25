
#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_double.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct Config2ndDoubleYZ
{
  using mparticles_t = PscMparticlesDouble;
  using mfields_t = PscMfieldsC;
};

#define CONFIG Config2ndDoubleYZ

#define DIM DIM_YZ
#define ORDER ORDER_2ND
#define PRTS PRTS_STAGGERED
//#define CACHE CACHE_EM_J
#define do_push_part do_push_part_2nd_yz
#define PROF_NAME "2nd_push_mprts_yz"
#include "push_part_common.c"

