
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#define DIM DIM_YZ

#include "inc_defs.h"

#define ORDER ORDER_1ST
#define CALC_J CALC_J_1VB_2D
#define EXT_PREPARE_SORT

#include "1vb.c"

template<typename dim_t>
using push_p_ops_1vb2_single2 = push_p_ops<push_p_config<mfields_single_t, dim_t>>;

// ======================================================================
// psc_push_particles: subclass "1vb2_single2"

struct psc_push_particles_ops_1vb2_single2 : psc_push_particles_ops {
  psc_push_particles_ops_1vb2_single2() {
    name                  = "1vb2_single2";
    push_mprts_yz         = push_p_ops_1vb2_single2<dim_yz>::push_mprts;
    particles_type        = PARTICLE_TYPE;
  }
} psc_push_particles_1vb2_single_ops;

