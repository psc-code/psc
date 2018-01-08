
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#define push_p_ops push_p_ops_1vbec_single_xz
#include "1vb/psc_push_particles_1vb.h"

#define DIM DIM_XYZ

#include "../inc_defs.h"

#define PUSH_DIM DIM_XZ
#define EM_CACHE_DIM DIM_XZ
#define CURR_CACHE_DIM DIM_XZ
#define ORDER ORDER_1ST
#define IP_VARIANT IP_VARIANT_EC
#define CALC_J CALC_J_1VB_SPLIT

#include "1vb.c"

