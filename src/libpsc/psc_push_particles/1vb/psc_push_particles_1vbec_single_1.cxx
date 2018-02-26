
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#define DIM DIM_1

#include "../inc_defs.h"
#include "../push_config.hxx"

#define EM_CACHE_DIM DIM_1
#define CURR_CACHE_DIM DIM_1
#define CALC_J CALC_J_1VB_VAR1

using push_p_conf = push_p_config<mparticles_t, mfields_t, dim_1, opt_ip_1st_ec, opt_order_1st,
				  opt_calcj_1vb_var1>;

#include "../1vb.c"

