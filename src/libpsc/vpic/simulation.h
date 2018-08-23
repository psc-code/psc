
#ifndef SIMULATION_H
#define SIMULATION_H

#include "psc_vpic_bits.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class VpicSimulation

template<class Mparticles, class MfieldsState, class MfieldsInterpolator, class MfieldsHydro,
	 class RP, class SimulationMixin, class DiagMixin>
struct VpicSimulation : SimulationMixin, DiagMixin
{
  using Particles = typename Mparticles::Particles;
  using Grid = typename Particles::Grid;
  using ParticleBcList = typename Particles::ParticleBcList;
  using MaterialList = typename MfieldsState::MaterialList;

  VpicSimulation()
    : SimulationMixin(),
      DiagMixin()
  {}

  // ----------------------------------------------------------------------
  // DiagMixin
  
  void newDiag(int interval)
  {
    DiagMixin::diagnostics_init(interval);
  }

  void setupDiag()
  {
    DiagMixin::diagnostics_setup();
  }

  void runDiag(Mparticles& mprts, MfieldsState& mflds, MfieldsInterpolator& interpolator,
	       MfieldsHydro& hydro, Int3 np)
  {
    DiagMixin::diagnostics_run(mprts, mflds, interpolator, hydro, np);
  }

  // ======================================================================
  // substeps of a time integration step

  // ----------------------------------------------------------------------
  // field_injection
  
  using SimulationMixin::field_injection;

  //private:
  Grid* vgrid_;
};

#endif



