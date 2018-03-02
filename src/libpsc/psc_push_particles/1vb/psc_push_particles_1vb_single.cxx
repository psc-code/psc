
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#define CALC_J CALC_J_1VB_2D

#include "../inc_defs.h"
#include "../push_config.hxx"

using push_p_conf = push_p_config<MparticlesSingle, MfieldsSingle,
				  InterpolateEM1st<Fields3d<MfieldsSingle::fields_t>, dim_yz>,
				  dim_yz, opt_order_1st,
				  Current1vb2d>;

#include "../1vb.c"

template<typename dim_t>
using push_p_ops_1vb_single = push_p_ops<push_p_conf, PushParticles1vb>;

using PushParticles_t = PushParticles_<push_p_ops_1vb_single>;
using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles_t>;
  
// ======================================================================
// psc_push_particles: subclass "1vb_single"

struct psc_push_particles_ops_1vb_single : psc_push_particles_ops {
  psc_push_particles_ops_1vb_single() {
    name                  = "1vb_single";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    particles_type        = PARTICLE_TYPE;
  }
} psc_push_particles_1vb_single_ops;

