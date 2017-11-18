
#ifndef VPIC_PUSH_PARTICLES_H
#define VPIC_PUSH_PARTICLES_H

#include "vpic_iface.h"

#include "simulation.h"

#include <vpic.h>

// ======================================================================
// vpic_push_particles

struct vpic_push_particles {
  InterpolatorArray* interpolator_array;
  AccumulatorArray* accumulator_array;
  int num_comm_round;

  void clear_accumulator_array();
  void advance_p(Particles *vmprts);
  void reduce_accumulator_array();
  void boundary_p(Particles *vmprts, FieldArray *vmflds);
  void unload_accumulator_array(FieldArray *vmflds);
  void load_interpolator_array(FieldArray *vmflds);
  void uncenter_p(Particles *vmprts);
};


#endif

