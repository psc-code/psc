
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM DIM_XZ
#define ORDER ORDER_1ST
#define PRTS PRTS_STAGGERED
#define IP_VARIANT IP_VARIANT_SFF
#include "push_part_common.c"

