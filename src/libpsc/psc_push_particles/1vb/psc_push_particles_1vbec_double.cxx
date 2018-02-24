
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"
#include "psc_push_particles_1vb.h"
#include "push_particles.hxx"

// ======================================================================
// psc_push_particles: subclass "1vbec_double"

template<typename dim_t>
using push_p_ops_1vbec_double = push_p_ops<push_p_config<mfields_c_t, dim_t>>;

using PushParticles_t = PushParticles_<push_p_ops_1vbec_double>;
  
struct psc_push_particles_ops_1vbec_double : psc_push_particles_ops {
  psc_push_particles_ops_1vbec_double() {
    name                  = "1vbec_double";
    size                  = sizeof(PushParticles_t);
    setup                 = PushParticles_t::setup;
    destroy               = PushParticles_t::destroy;
    push_mprts_xyz        = PushParticles_t::push_mprts_xyz;
    push_mprts_yz         = PushParticles_t::push_mprts_yz;
    push_mprts_1          = PushParticles_t::push_mprts_1;
    //stagger_mprts_1      = push_p_ops_1vbec_double<dim_1>::stagger_mprts;
    particles_type        = PARTICLE_TYPE;
  }
} psc_push_particles_1vbec_double_ops;

