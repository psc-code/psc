
#ifndef SIMULATION_H
#define SIMULATION_H

#include "hydro_array.h"
#include "material.h"
#include "rng.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class VpicSimulation

template<class FieldArrayOps, class ParticlesOps, class InterpolatorOps,
	 class AccumulatorOps, class SimulationBase>
struct VpicSimulation : FieldArrayOps, ParticlesOps, InterpolatorOps, AccumulatorOps
{
  typedef typename FieldArrayOps::FieldArray FieldArray;
  typedef typename ParticlesOps::Particles Particles;
  typedef typename InterpolatorOps::Interpolator Interpolator;
  typedef typename AccumulatorOps::Accumulator Accumulator;
  
  VpicSimulation(SimulationBase *sim_base)
    : ParticlesOps(reinterpret_cast<vpic_simulation*>(sim_base)),
      num_comm_round_(sim_base->num_comm_round),
      grid_(sim_base->getGrid()),
      material_list_(sim_base->getMaterialList()),
      field_array_(sim_base->getFieldArray()),
      interpolator_(sim_base->getInterpolator()),
      accumulator_(sim_base->getAccumulator()),
      hydro_array_(sim_base->hydro_array),
      particles_(sim_base->getParticles()),
      sim_base_(sim_base)
  {
  }

  void set_params(int num_step, int status_interval,
		  int sync_shared_interval, int clean_div_e_interval,
		  int clean_div_b_interval)
  {
    sim_base_->setParams(num_step, status_interval, sync_shared_interval,
			 clean_div_e_interval, clean_div_b_interval);
  }
  
  void setup_grid(double dx[3], double dt, double cvac, double eps0)
  {
    grid_.setup(dx, dt, cvac, eps0);
  }

  void define_periodic_grid(double xl[3], double xh[3], int gdims[3], int np[3])
  {
    sim_base_->setTopology(np[0], np[1], np[2]);
    grid_.partition_periodic_box(xl, xh, gdims, np);
  }

  void set_domain_field_bc(int boundary, int bc)
  {
    int fbc;
    switch (bc) {
    case BND_FLD_CONDUCTING_WALL: fbc = pec_fields   ; break;
    case BND_FLD_ABSORBING:       fbc = absorb_fields; break;
    default: assert(0);
    }
    grid_.set_fbc(boundary, fbc);
  }

  void set_domain_particle_bc(int boundary, int bc)
  {
    int pbc;
    switch (bc) {
    case BND_PART_REFLECTING: pbc = reflect_particles; break;
    case BND_PART_ABSORBING:  pbc = absorb_particles ; break;
    default: assert(0);
    }
    grid_.set_pbc(boundary, pbc);
  }

  struct material *define_material(const char *name,
						   double eps, double mu=1.,
						   double sigma=0., double zeta=0.)
  {
    return material_list_.append(material(name,
					  eps,   eps,   eps,
					  mu,    mu,    mu,
					  sigma, sigma, sigma,
					  zeta,  zeta,  zeta));
  }

  void define_field_array(double damp)
  {
    grid_t *g = grid_.getGrid_t();
 
    assert(g->nx && g->ny && g->ny);
    assert(!material_list_.empty());
  
    field_array_ = new FieldArray(grid_, material_list_, damp);
    interpolator_ = new Interpolator(g);
    accumulator_ = new Accumulator(g);
    hydro_array_ = new HydroArray(grid_);
 
    // Pre-size communications buffers. This is done to get most memory
    // allocation over with before the simulation starts running
  
    int nx1 = g->nx+1, ny1 = g->ny+1, nz1 = g->nz+1;
    grid_.mp_size_recv_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_.mp_size_recv_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_.mp_size_recv_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
    grid_.mp_size_recv_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
    grid_.mp_size_recv_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
    grid_.mp_size_recv_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
  
    grid_.mp_size_send_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_.mp_size_send_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
    grid_.mp_size_send_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
    grid_.mp_size_send_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
    grid_.mp_size_send_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
    grid_.mp_size_send_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
  }

  species_t* define_species(const char *name, double q, double m,
			    double max_local_np, double max_local_nm,
			    double sort_interval, double sort_out_of_place)
  {
    // Compute a reasonble number of movers if user did not specify
    // Based on the twice the number of particles expected to hit the boundary
    // of a wpdt=0.2 / dx=lambda species in a 3x3x3 domain
    if (max_local_nm < 0) {
      max_local_nm = 2 * max_local_np / 25;
      if (max_local_nm<16*(MAX_PIPELINE+1))
	max_local_nm = 16*(MAX_PIPELINE+1);
    }
    return particles_.append(species(name, (float)q, (float)m,
				     (int)max_local_np, (int)max_local_nm,
				     (int)sort_interval, (int)sort_out_of_place,
				     grid_.getGrid_t()));
  }

  void collision_run()
  {
    sim_base_->collision_run();
  }
  
  void emitter()
  {
    sim_base_->emitter();
  }

  void current_injection()
  {
    sim_base_->user_current_injection();
  }

  void field_injection()
  {
    sim_base_->user_field_injection();
  }

  void moments_run(HydroArray *hydro_array, Particles *vmprts, int kind)
  {
    // This relies on load_interpolator_array() having been called earlier
    clear_hydro_array(hydro_array);

    typename Particles::Iter sp = vmprts->find_id(kind);
    accumulate_hydro_p(hydro_array, &*sp, interpolator_);
    
    synchronize_hydro_array(hydro_array);
  }

  void newDiag(int interval)
  {
    sim_base_->diagInit(interval);
  }

  void setupDiag()
  {
    sim_base_->diagSetup();
  }

  void runDiag()
  {
    sim_base_->diagRun();
  }
    
  int num_comm_round_;
  
  RngPool rng_pool;

  //private:
  Grid grid_;
  MaterialList& material_list_;
  FieldArray*& field_array_;
  Interpolator*& interpolator_;
  Accumulator*& accumulator_;
  hydro_array_t*& hydro_array_;
  Particles& particles_;

  SimulationBase *sim_base_;
};

#endif



