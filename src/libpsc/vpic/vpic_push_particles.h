
#ifndef VPIC_PUSH_PARTICLES_H
#define VPIC_PUSH_PARTICLES_H

#include "vpic_iface.h"

#include <vpic.h>

// ======================================================================
// vpic_push_particles

struct vpic_push_particles {
  interpolator_array_t *interpolator_array;
  accumulator_array_t *accumulator_array;
  int num_comm_round;

  void clear_accumulator_array();
  void advance_p(vpic_mparticles *vmprts);
  void reduce_accumulator_array();
  void boundary_p(vpic_mparticles *vmprts, FieldArray *vmflds);
  void unload_accumulator_array(FieldArray *vmflds);
  void load_interpolator_array(FieldArray *vmflds);
  void uncenter_p(vpic_mparticles *vmprts);
};


#endif

