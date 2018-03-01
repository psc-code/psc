
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#define DIM DIM_YZ

#include "../inc_defs.h"
#include "../push_config.hxx"

using push_p_conf = Config1vbecSingle<dim_yz>;

#include "../1vb.c"

