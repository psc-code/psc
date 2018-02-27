
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#define DIM DIM_YZ

#include "../inc_defs.h"
#include "../push_config.hxx"

#define CALC_J CALC_J_1VB_VAR1

using push_p_conf = push_p_config<MparticlesDouble, MfieldsC, dim_yz, opt_ip_1st_ec, opt_order_1st,
				  opt_calcj_1vb_var1>;

#include "../1vb.c"

