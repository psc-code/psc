
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#define DIM DIM_XYZ

#include "../inc_defs.h"
#include "../push_config.hxx"

#define PUSH_DIM DIM_XZ
#define EM_CACHE_DIM DIM_XZ
#define CURR_CACHE_DIM DIM_XZ
#define CALC_J CALC_J_1VB_SPLIT

using push_p_conf = Config1vbecSingleXZ;

#include "../1vb.c"

