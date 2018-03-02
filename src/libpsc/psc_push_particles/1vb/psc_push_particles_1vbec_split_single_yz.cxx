
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_push_particles_1vb.h"

#include "../inc_defs.h"
#include "../push_config.hxx"

#define CALC_J CALC_J_1VB_SPLIT

#include "../1vb.c"

template struct push_p_ops<push_p_config<MparticlesSingle, MfieldsSingle,
					 InterpolateEM1vbec<Fields3d<MfieldsSingle::fields_t>, dim_yz>,
					 dim_yz, opt_order_1st,
					 Current1vbSplit>,
			   PushParticles1vb>;

