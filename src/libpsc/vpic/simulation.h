
#ifndef SIMULATION_H
#define SIMULATION_H

#include "psc_vpic_bits.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class VpicSimulation

template<class P, class RP, class SimulationMixin, class DiagMixin>
struct VpicSimulation : SimulationMixin, DiagMixin
{
  typedef P Particles;
  typedef typename Particles::Grid Grid;
  typedef typename Particles::FieldArray FieldArray;
  typedef typename Particles::Interpolator Interpolator;
  typedef typename Particles::Accumulator Accumulator;
  typedef typename Particles::HydroArray HydroArray;
  typedef typename Particles::Species Species;
  typedef typename Particles::ParticleBcList ParticleBcList;
  typedef typename FieldArray::MaterialList MaterialList;
  typedef typename MaterialList::Material Material;
  typedef RP RngPool;

  VpicSimulation()
    : SimulationMixin(),
      DiagMixin(),
      num_comm_round_(3),
      grid_(SimulationMixin::getGrid()),
      particle_bc_list_(SimulationMixin::getParticleBcList())
  {
  }

  void define_periodic_grid(double xl[3], double xh[3], const int gdims[3], const int np[3])
  {
    np_[0] = np[0]; np_[1] = np[1]; np_[2] = np[2];
    SimulationMixin::setTopology(np[0], np[1], np[2]);
    grid_->partition_periodic_box(xl, xh, gdims, np);
  }

  void set_domain_field_bc(int boundary, int bc)
  {
    int fbc;
    switch (bc) {
    case BND_FLD_CONDUCTING_WALL: fbc = Grid::pec_fields   ; break;
    case BND_FLD_ABSORBING:       fbc = Grid::absorb_fields; break;
    default: assert(0);
    }
    grid_->set_fbc(boundary, fbc);
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

  Material *define_material(const char *name,
			    double eps, double mu=1.,
			    double sigma=0., double zeta=0.)
  {
    Material *m = material_list_.create(name,
					eps,   eps,   eps,
					mu,    mu,    mu,
					sigma, sigma, sigma,
					zeta,  zeta,  zeta);
    return material_list_.append(m);
  }

  void grid_setup_communication()
  {
    assert(grid_->nx && grid_->ny && grid_->ny);
  
    // Pre-size communications buffers. This is done to get most memory
    // allocation over with before the simulation starts running
    // FIXME, this isn't a great place. First, we shouldn't call mp
    // functions (semi-)directly. 2nd, whether we need these buffers depends
    // on b.c., which aren't yet known.

    // FIXME, this really isn't a good place to do this, as it requires layer breaking knowledge of
    // which communication will need the largest buffers...
    int nx1 = grid_->nx+1, ny1 = grid_->ny+1, nz1 = grid_->nz+1;
    grid_->mp_size_recv_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(typename HydroArray::Element));
    grid_->mp_size_recv_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(typename HydroArray::Element));
    grid_->mp_size_recv_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(typename HydroArray::Element));
    grid_->mp_size_recv_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(typename HydroArray::Element));
    grid_->mp_size_recv_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(typename HydroArray::Element));
    grid_->mp_size_recv_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(typename HydroArray::Element));
  
    grid_->mp_size_send_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(typename HydroArray::Element));
    grid_->mp_size_send_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(typename HydroArray::Element));
    grid_->mp_size_send_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(typename HydroArray::Element));
    grid_->mp_size_send_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(typename HydroArray::Element));
    grid_->mp_size_send_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(typename HydroArray::Element));
    grid_->mp_size_send_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(typename HydroArray::Element));
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

  void runDiag(Particles& particles, FieldArray& fa, Interpolator& ia, HydroArray& ha)
  {
    DiagMixin::diagnostics_run(fa, particles, ia, ha, np_);
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
  ParticleBcList& particle_bc_list_;

  int np_[3];
};

#endif



