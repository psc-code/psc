
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#define DIM DIM_1

#include "../inc_defs.h"

#define EM_CACHE_DIM DIM_1
#define CURR_CACHE_DIM DIM_1
#define ORDER ORDER_1ST
#define IP_VARIANT IP_VARIANT_EC
#define CALC_J CALC_J_1VB_VAR1

#include "../1vb.c"

