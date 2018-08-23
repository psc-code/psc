
#ifndef SIMULATION_H
#define SIMULATION_H

#include "psc_vpic_bits.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class VpicSimulation

template<class Mparticles, class MfieldsState, class MfieldsInterpolator, class MfieldsHydro,
	 class RP, class DiagMixin>
struct VpicSimulation : DiagMixin
{
  using Grid = typename Mparticles::Grid;

  //private:
  Grid* vgrid_;
};

#endif



