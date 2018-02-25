
#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_double.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "push_config.hxx"

#define CONFIG Config2ndDoubleYZ

#define DIM DIM_YZ
#define ORDER ORDER_2ND
#define PRTS PRTS_STAGGERED
#define CACHE CACHE_EM_J
#include "push_part_common.c"

