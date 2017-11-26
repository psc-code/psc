
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
  VpicInterpolatorBase* interpolator;
  VpicAccumulatorBase* accumulator;
  int num_comm_round;

  void push_mprts(Particles *vmprts, FieldArray *vmflds);
  void prep(Particles *vmprts, FieldArray *vmflds);
};

#endif

