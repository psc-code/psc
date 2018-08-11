
#ifndef SIMULATION_H
#define SIMULATION_H

#include "psc_vpic_bits.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class VpicSimulation

template<class P, class FieldArray, class Interpolator, class HydroArray, class RP, class SimulationMixin, class DiagMixin>
struct VpicSimulation : SimulationMixin, DiagMixin
{
  typedef P Particles;
  typedef typename Particles::Grid Grid;
  typedef typename Particles::Species Species;
  typedef typename Particles::ParticleBcList ParticleBcList;
  typedef typename FieldArray::MaterialList MaterialList;
  typedef typename MaterialList::Material Material;
  typedef RP RngPool;

  VpicSimulation()
    : SimulationMixin(),
      DiagMixin(),
      num_comm_round_(3),
      grid_(SimulationMixin::getGrid())
  {
  }

  void set_domain_particle_bc(int boundary, int bc)
  {
    int pbc;
    switch (bc) {
    case BND_PRT_REFLECTING: pbc = Grid::reflect_particles; break;
    case BND_PRT_ABSORBING:  pbc = Grid::absorb_particles ; break;
    default: assert(0);
    }
    grid_->set_pbc(boundary, pbc);
  }

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

  void runDiag(Particles& particles, FieldArray& fa, Interpolator& ia, HydroArray& ha, Int3 np)
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

  int num_comm_round_;
  
  //private:
  Grid*& grid_;
  MaterialList material_list_;
  ParticleBcList particle_bc_list_;
};

#endif



