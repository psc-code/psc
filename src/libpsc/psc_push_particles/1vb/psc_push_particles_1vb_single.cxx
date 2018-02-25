
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#define DIM DIM_YZ

#include "../inc_defs.h"

#define ORDER ORDER_1ST
#define CALC_J CALC_J_1VB_2D

#include "../1vb.c"

template<typename dim_t>
using push_p_ops_1vb_single = push_p_ops<push_p_config<mparticles_single_t, mfields_single_t, dim_t, opt_order_1st, opt_calcj_1vb_2d>>;

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

