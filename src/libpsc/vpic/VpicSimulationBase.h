
#ifndef VPIC_SIMULATION_BASE_H
#define VPIC_SIMULATION_BASE_H

#include "VpicParticlesBase.h"
#include "VpicDiag.h"

// FIXME, the casts below will happily do the wrong thing if
// this the underlying base types aren't vpic / vpic-compatible layout

template<class Particles>
class VpicSimulationMixin : protected vpic_simulation
{
  typedef typename Particles::Grid Grid;

public:
  VpicSimulationMixin()
  {
    extern vpic_simulation *simulation;
    assert(!simulation);
    simulation = this;
  }
};


#endif
