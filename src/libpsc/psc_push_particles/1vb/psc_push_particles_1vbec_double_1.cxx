
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#define CALC_J CALC_J_1VB_VAR1

#include "../inc_defs.h"
#include "../push_config.hxx"

#include "../1vb.c"

template struct push_p_ops<Config1vbecDouble<dim_1>>;
