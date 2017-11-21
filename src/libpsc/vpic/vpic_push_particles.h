
#ifndef VPIC_PUSH_PARTICLES_H
#define VPIC_PUSH_PARTICLES_H

#include "vpic_iface.h"

#include "simulation.h"

#include <vpic.h>

// ======================================================================
// vpic_push_particles

struct vpic_push_particles {
  vpic_push_particles(Simulation *sim);
      
  Simulation *sim_;
  VpicInterpolator* interpolator;
  VpicAccumulator* accumulator_array;
  int num_comm_round;

  void stagger_mprts(Particles *vmprts, FieldArray *vmflds);
  void push_mprts(Particles *vmprts, FieldArray *vmflds);
  void prep(Particles *vmprts, FieldArray *vmflds);
};

#endif

