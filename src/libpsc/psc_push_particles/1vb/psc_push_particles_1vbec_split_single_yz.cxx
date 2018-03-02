
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#include "../inc_defs.h"
#include "../push_config.hxx"

#define CALC_J CALC_J_1VB_SPLIT

using push_p_conf = push_p_config<MparticlesSingle, MfieldsSingle,
				  InterpolateEM<Fields3d<MfieldsSingle::fields_t>, opt_ip_1st_ec, dim_yz>,
				  dim_yz, opt_order_1st,
				  Current1vbSplit>;

#include "../1vb.c"

