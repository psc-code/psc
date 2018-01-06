
#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_double.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM DIM_YZ
#define ORDER ORDER_2ND
#define PRTS PRTS_STAGGERED
//#define CACHE CACHE_EM_J
#define psc_push_particles_push_mprts psc_push_particles_2nd_double_push_mprts_yz
#define do_push_part do_push_part_2nd_yz
#define PROF_NAME "2nd_push_mprts_yz"
#include "push_part_common.c"

