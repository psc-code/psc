
#include "psc_1p5_c.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM DIM_YZ
#define ORDER ORDER_1P5
#define psc_push_particles_push_mprts psc_push_particles_1p5_c_push_mprts_yz
#define do_push_part do_push_part_1p5_yz
#define PROF_NAME "push_mprts_1p5_yz"
#include "push_part_common.c"

