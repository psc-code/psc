
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#define DIM DIM_YZ

#include "inc_defs.h"
#include "push_config.hxx"

#define ORDER ORDER_1ST
#define CALC_J CALC_J_1VB_2D

using push_p_conf = push_p_config<PscMparticlesSingle, PscMfieldsSingle, dim_yz,
				  opt_ip_1st, opt_order_1st, opt_calcj_1vb_2d, opt_ext_prepare_sort>;

#include "1vb.c"

template<typename dim_t>
using push_p_ops_1vb2_single2 = push_p_ops<push_p_conf>;

using PushParticles_t = PushParticles_<push_p_ops_1vb2_single2>;

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles_t>;

// ======================================================================
// psc_push_particles: subclass "1vb2_single2"

struct psc_push_particles_ops_1vb2_single2 : psc_push_particles_ops {
  psc_push_particles_ops_1vb2_single2() {
    name                  = "1vb2_single2";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    particles_type        = PARTICLE_TYPE;
  }
} psc_push_particles_1vb2_single_ops;

