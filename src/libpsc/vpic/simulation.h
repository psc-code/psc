
#ifndef SIMULATION_H
#define SIMULATION_H

#include "psc_vpic_bits.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class VpicSimulation

template<class Particles, class FieldArray, class Interpolator, class MfieldsHydro,
	 class RP, class SimulationMixin, class DiagMixin>
struct VpicSimulation : SimulationMixin, DiagMixin
{
  using Grid = typename Particles::Grid;
  using ParticleBcList = typename Particles::ParticleBcList;
  using MaterialList = typename FieldArray::MaterialList;

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

  void runDiag(Particles& particles, FieldArray& fa, Interpolator& ia, MfieldsHydro& ha, Int3 np)
  {
    DiagMixin::diagnostics_run(fa, particles, ia, ha, np);
  }

  // ======================================================================
  // substeps of a time integration step

  // ----------------------------------------------------------------------
  // collision_run
  
  using SimulationMixin::collision_run;

  // ----------------------------------------------------------------------
  // field_injection
  
  using SimulationMixin::field_injection;

  //private:
  Grid* vgrid_;
};

#endif



