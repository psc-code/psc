
#include "psc_push_particles_private.h"

#include "psc_particles_as_single_by_block.h"
#include "psc_fields_as_single.h"
#include "1vb/psc_push_particles_1vb.h"

// ======================================================================
// psc_push_particles: subclass "1vbec_single_by_block"

template<typename dim_t>
using push_p_ops_1vb_single_by_block = push_p_ops<push_p_config<mfields_single_t, dim_t>>;

struct psc_push_particles_ops psc_push_particles_1vbec_single_by_block_ops = {
  .name                  = "1vbec_single_by_block",
  .push_mprts_yz         = push_p_ops_1vb_single_by_block<dim_yz>::push_mprts,
  .push_mprts_xyz        = push_p_ops_1vb_single_by_block<dim_xyz>::push_mprts,
  .particles_type        = PARTICLE_TYPE,
};

