
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#define DIM DIM_YZ

#include "../inc_defs.h"
#include "../push_config.hxx"

#define CALC_J CALC_J_1VB_SPLIT

using push_p_conf = push_p_config<MparticlesSingle, MfieldsSingle, dim_yz, opt_ip_1st_ec, opt_order_1st,
				  Current1vbSplit, opt_calcj_1vb_var1>;

#include "../1vb.c"

