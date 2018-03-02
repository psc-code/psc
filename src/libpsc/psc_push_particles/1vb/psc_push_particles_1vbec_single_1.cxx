
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#include "../inc_defs.h"
#include "../push_config.hxx"

#define EM_CACHE_DIM DIM_1
#define CURR_CACHE_DIM DIM_1
#define CALC_J CALC_J_1VB_VAR1

using push_p_conf = Config1vbecSingle1;

#include "../1vb.c"

