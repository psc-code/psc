
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_single.h"

#define push_p_ops push_p_ops_1vbec_single
#include "psc_push_particles_1vb.h"

// ======================================================================
// psc_push_particles: subclass "1vbec_single"

template<typename dim_t>
using push_p_ops_1vbec_single_ = push_p_ops_1vbec_single<push_p_config<mfields_single_t, dim_t>>;

// FIXME, special hack... for xyz_xz
template<typename C>
struct push_p_ops_1vbec_single_xz
{
  static void push_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base);
  static void stagger_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base);
};

struct psc_push_particles_ops_1vbec_single : psc_push_particles_ops {
  psc_push_particles_ops_1vbec_single() {
    name                  = "1vbec_single";
    push_mprts_xyz        = push_p_ops_1vbec_single_<dim_xyz>::push_mprts;
    push_mprts_xz         = push_p_ops_1vbec_single_xz<push_p_config<mfields_single_t, dim_xyz>>::push_mprts;
    push_mprts_yz         = push_p_ops_1vbec_single_<dim_yz>::push_mprts;
    push_mprts_1          = push_p_ops_1vbec_single_<dim_1>::push_mprts;
    stagger_mprts_yz      = push_p_ops_1vbec_single_<dim_yz>::stagger_mprts;
    particles_type        = PARTICLE_TYPE;
  }
} psc_push_particles_1vbec_single_ops;

