
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct Config1stXZ;

#define CONFIG Config1stXZ

#define DIM DIM_XZ
#define ORDER ORDER_1ST
#define PRTS PRTS_STAGGERED
#include "push_part_common.c"

