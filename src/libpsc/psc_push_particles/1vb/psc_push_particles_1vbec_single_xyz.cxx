
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#define DIM DIM_XYZ
#define CALC_J CALC_J_1VB_VAR1

#include "../inc_defs.h"
#include "../push_config.hxx"

using push_p_conf = Config1vbecSingle<dim_1>;

#include "../1vb.c"

