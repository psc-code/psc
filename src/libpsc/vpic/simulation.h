
#ifndef SIMULATION_H
#define SIMULATION_H

#include "material.h"
#include "rng.h"

#include "species_advance/species_advance.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class VpicSimulation

template<class P, class ParticlesOps, class RP, class SimulationMixin, class DiagMixin>
struct VpicSimulation : SimulationMixin, ParticlesOps, DiagMixin
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

  using SimulationMixin::collision_run;
  using SimulationMixin::field_injection;

  VpicSimulation()
    : SimulationMixin(),
      ParticlesOps(),
      DiagMixin(),
      num_comm_round_(3),
      grid_(SimulationMixin::getGrid()),
      material_list_(SimulationMixin::getMaterialList()),
      field_array_(SimulationMixin::getFieldArray()),
      interpolator_(SimulationMixin::getInterpolator()),
      accumulator_(SimulationMixin::getAccumulator()),
      hydro_array_(SimulationMixin::getHydroArray()),
      particles_(SimulationMixin::getParticles()),
      particle_bc_list_(SimulationMixin::getParticleBcList())
  {
  }

  void setup_grid(double dx[3], double dt, double cvac, double eps0)
  {
    grid_->setup(dx, dt, cvac, eps0);
  }

  void define_periodic_grid(double xl[3], double xh[3], int gdims[3], int np[3])
  {
    np_[0] = np[0]; np_[1] = np[1]; np_[2] = np[2];
    SimulationMixin::setTopology(np[0], np[1], np[2]);
    grid_->partition_periodic_box(xl, xh, gdims, np);
  }

  void set_domain_field_bc(int boundary, int bc)
  {
    int fbc;
    switch (bc) {
    case BND_FLD_CONDUCTING_WALL: fbc = pec_fields   ; break;
    case BND_FLD_ABSORBING:       fbc = absorb_fields; break;
    default: assert(0);
    }
    grid_->set_fbc(boundary, fbc);
  }

  void set_domain_particle_bc(int boundary, int bc)
  {
    int pbc;
    switch (bc) {
    case BND_PART_REFLECTING: pbc = reflect_particles; break;
    case BND_PART_ABSORBING:  pbc = absorb_particles ; break;
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

  void define_field_array(double damp)
  {
    assert(grid_->nx && grid_->ny && grid_->ny);
    assert(!material_list_.empty());
  
    field_array_ = FieldArray::create(grid_, material_list_, damp);
    interpolator_ = Interpolator::create(grid_);
    accumulator_ = Accumulator::create(grid_);
    hydro_array_ = HydroArray::create(grid_);
 
    // Pre-size communications buffers. This is done to get most memory
    // allocation over with before the simulation starts running
    // FIXME, this isn't a great place. First, we shouldn't call mp
    // functions (semi-)directly. 2nd, whether we need these buffers depends
    // on b.c., which aren't yet known.
  
    int nx1 = grid_->nx+1, ny1 = grid_->ny+1, nz1 = grid_->nz+1;
    grid_->mp_size_recv_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_->mp_size_recv_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_->mp_size_recv_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
    grid_->mp_size_recv_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
    grid_->mp_size_recv_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
    grid_->mp_size_recv_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
  
    grid_->mp_size_send_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_->mp_size_send_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_->mp_size_send_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
    grid_->mp_size_send_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
    grid_->mp_size_send_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
    grid_->mp_size_send_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
  }

  Species* define_species(const char *name, double q, double m,
			  double max_local_np, double max_local_nm,
			  double sort_interval, double sort_out_of_place)
  {
    // Compute a reasonble number of movers if user did not specify
    // Based on the twice the number of particles expected to hit the boundary
    // of a wpdt=0.2 / dx=lambda species in a 3x3x3 domain
    if (max_local_nm < 0) {
      max_local_nm = 2 * max_local_np / 25;
      if (max_local_nm < 16*(MAX_PIPELINE+1))
	max_local_nm = 16*(MAX_PIPELINE+1);
    }
    Species *sp = particles_.create(name, q, m, max_local_np, max_local_nm,
				    sort_interval, sort_out_of_place, grid_);
    return particles_.append(sp);
  }

  void moments_run(HydroArray *hydro_array, Particles *vmprts, int kind)
  {
    // This relies on load_interpolator_array() having been called earlier
    hydro_array->clear();

    typename Particles::const_iterator sp = vmprts->find(kind);
    vmprts->accumulate_hydro_p(*hydro_array, sp, *interpolator_);
    
    hydro_array->synchronize();
  }

  void uncenter_p(Particles *vmprts, FieldArray *vmflds)
  {
    if (!vmprts->empty()) {
      interpolator_->load(*vmflds);
      
      for (auto sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
	TIC vmprts->uncenter_p(&*sp, *interpolator_); TOC(uncenter_p, 1);
      }
    }
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

  void runDiag()
  {
    DiagMixin::diagnostics_run(*field_array_, particles_, *interpolator_, *hydro_array_, np_);
  }

  // ----------------------------------------------------------------------
  // push_mprts

  void push_mprts(Particles& vmprts, FieldArray& vmflds)
  {
    // For this to work, interpolator_ needs to have been set from vmflds E/B before,
    // ie., we're not using vmflds for E and B here at all.
    
    // At this point, fields are at E_0 and B_0 and the particle positions
    // are at r_0 and u_{-1/2}.  Further the mover lists for the particles should
    // empty and all particles should be inside the local computational domain.
    // Advance the particle lists.
    if (!vmprts.empty()) {
      TIC accumulator_->clear(); TOC(clear_accumulators, 1);
      this->advance_p(vmprts, *accumulator_, *interpolator_);
    }

    // Because the partial position push when injecting aged particles might
    // place those particles onto the guard list (boundary interaction) and
    // because advance_p requires an empty guard list, particle injection must
    // be done after advance_p and before guard list processing. Note:
    // user_particle_injection should be a stub if sl_ is empty.
    this->emitter();

    // This should be after the emission and injection to allow for the
    // possibility of thread parallelizing these operations
    if (!vmprts.empty()) {
      TIC accumulator_->reduce(); TOC(reduce_accumulators, 1);
    }
    
    // At this point, most particle positions are at r_1 and u_{1/2}. Particles
    // that had boundary interactions are now on the guard list. Process the
    // guard lists. Particles that absorbed are added to rhob (using a corrected
    // local accumulation).
    TIC
      for(int round = 0; round < num_comm_round_; round++) {
	this->boundary_p(particle_bc_list_, vmprts, vmflds, *accumulator_);
      } TOC(boundary_p, num_comm_round_);

    for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
      if (sp->nm) {
	LOG_WARN("Removing %i particles associated with unprocessed %s movers (increase num_comm_round)",
		 sp->nm, sp->name);
      }
      // Drop the particles that have unprocessed movers due to a user defined
      // boundary condition. Particles of this type with unprocessed movers are
      // in the list of particles and move_p has set the voxel in the particle to
      // 8*voxel + face. This is an incorrect voxel index and in many cases can
      // in fact go out of bounds of the voxel indexing space. Removal is in
      // reverse order for back filling. Particle charge is accumulated to the
      // mesh before removing the particle.
      int nm = sp->nm;
      particle_mover_t * RESTRICT ALIGNED(16)  pm = sp->pm + sp->nm - 1;
      particle_t * RESTRICT ALIGNED(128) p0 = sp->p;
      for (; nm; nm--, pm--) {
	int i = pm->i; // particle index we are removing
	p0[i].i >>= 3; // shift particle voxel down
	// accumulate the particle's charge to the mesh
	this->accumulate_rhob(vmflds, p0 + i, sp->q);
	p0[i] = p0[sp->np - 1]; // put the last particle into position i
	sp->np--; // decrement the number of particles
      }
      sp->nm = 0;
    }

    // At this point, all particle positions are at r_1 and u_{1/2}, the
    // guard lists are empty and the accumulators on each processor are current.
    // Convert the accumulators into currents.
    TIC vmflds.clear_jf(); TOC(clear_jf, 1);
    if (!vmprts.empty()) {
      TIC accumulator_->unload(vmflds); TOC(unload_accumulator, 1);
    }
    TIC vmflds.synchronize_jf(); TOC(synchronize_jf, 1);

    // At this point, the particle currents are known at jf_{1/2}.
    // Let the user add their own current contributions. It is the users
    // responsibility to insure injected currents are consistent across domains.
    // It is also the users responsibility to update rhob according to
    // rhob_1 = rhob_0 + div juser_{1/2} (corrected local accumulation) if
    // the user wants electric field divergence cleaning to work.
    TIC this->current_injection(); TOC(user_current_injection, 1);
  }
  
  void push_mprts_prep(FieldArray& vmflds)
  {
    // At end of step:
    // Fields are updated ... load the interpolator for next time step and
    // particle diagnostics in user_diagnostics if there are any particle
    // species to worry about
    
    if (!particles_.empty()) {
      interpolator_->load(vmflds);
    }
  }

  int num_comm_round_;
  
  RngPool rng_pool;

  //private:
  Grid*& grid_;
  MaterialList& material_list_;
  FieldArray*& field_array_;
  Interpolator*& interpolator_;
  Accumulator*& accumulator_;
  HydroArray*& hydro_array_;
  Particles& particles_;
  ParticleBcList& particle_bc_list_;

  int np_[3];
};

#endif



