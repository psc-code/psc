
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#include "../inc_defs.h"
#include "../push_config.hxx"

#define PUSH_DIM DIM_XZ
#define EM_CACHE_DIM DIM_XZ
#define CURR_CACHE_DIM DIM_XZ
#define CALC_J CALC_J_1VB_SPLIT

#include "../1vb.c"

template struct push_p_ops<Config1vbecSingleXZ>;

