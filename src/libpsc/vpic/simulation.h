
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag

#include "field_array.h"
#include "rng.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class Simulation

struct Simulation {
  Simulation();
  ~Simulation();

  void setup_grid(double dx[3], double dt, double cvac, double eps0);
  void define_periodic_grid(double xl[3], double xh[3], int gdims[3],
			    int np[3]);
  void set_domain_field_bc(int boundary, int bc);
  void set_domain_particle_bc(int boundary, int bc);

  struct material *define_material(const char *name, double eps, double mu,
				   double sigma, double zeta);
  field_array_t *new_field_array(float damp=0.);
  void define_field_array(double damp);
  species_t* define_species(const char *name, double q, double m,
			    double max_local_np, double max_local_nm,
			    double sort_interval, double sort_out_of_place);
    
  RngPool rng_pool;

  //private:
  globals_diag *pDiag_;
  Grid grid_;
  MaterialList material_list_;
  field_array_t*& field_array_;
  interpolator_array_t*& interpolator_array_;
  accumulator_array_t*& accumulator_array_;
  hydro_array_t*& hydro_array_;

  species_t*& species_list_;
};

inline Simulation::Simulation()
  : grid_(simulation->grid),
    material_list_(simulation->material_list),
    field_array_(simulation->field_array),
    interpolator_array_(simulation->interpolator_array),
    accumulator_array_(simulation->accumulator_array),
    hydro_array_(simulation->hydro_array),
    species_list_(simulation->species_list)
{
}

inline Simulation::~Simulation()
{
  delete pDiag_;
}

inline void Simulation::setup_grid(double dx[3], double dt, double cvac, double eps0)
{
  grid_.setup(dx, dt, cvac, eps0);
}

inline void Simulation::define_periodic_grid(double xl[3], double xh[3], int gdims[3],
					     int np[3])
{
  simulation->px = size_t(np[0]);
  simulation->py = size_t(np[1]);
  simulation->pz = size_t(np[2]);
  grid_.partition_periodic_box(xl, xh, gdims, np);
}

inline void Simulation::set_domain_field_bc(int boundary, int bc)
{
  int fbc;
  switch (bc) {
  case BND_FLD_CONDUCTING_WALL: fbc = pec_fields   ; break;
  case BND_FLD_ABSORBING:       fbc = absorb_fields; break;
  default: assert(0);
  }
  grid_.set_fbc(boundary, fbc);
}

inline void Simulation::set_domain_particle_bc(int boundary, int bc)
{
  int pbc;
  switch (bc) {
  case BND_PART_REFLECTING: pbc = reflect_particles; break;
  case BND_PART_ABSORBING:  pbc = absorb_particles ; break;
  default: assert(0);
  }
  grid_.set_pbc(boundary, pbc);
}

inline struct material *Simulation::define_material(const char *name,
						    double eps, double mu=1.,
						    double sigma=0., double zeta=0.)
{
  return material_list_.append(material(name,
					eps,   eps,   eps,
					mu,    mu,    mu,
					sigma, sigma, sigma,
					zeta,  zeta,  zeta));
}

inline void Simulation::define_field_array(double damp)
{
  grid_t *grid = grid_.g_;
 
  assert(grid_.g_->nx && grid_.g_->ny && grid_.g_->ny);
  assert(!material_list_.empty());
  
  field_array_ = new FieldArray(grid_, material_list_, damp);
  interpolator_array_ = ::new_interpolator_array(grid);
  accumulator_array_ = ::new_accumulator_array(grid);
  hydro_array_ = ::new_hydro_array(grid);
 
  // Pre-size communications buffers. This is done to get most memory
  // allocation over with before the simulation starts running
  
  int nx1 = grid->nx+1, ny1 = grid->ny+1, nz1 = grid->nz+1;
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

inline species_t* Simulation::define_species(const char *name, double q, double m,
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
  return append_species(species(name, (float)q, (float)m,
				(int)max_local_np, (int)max_local_nm,
				(int)sort_interval, (int)sort_out_of_place,
				grid_.g_), &species_list_);
}


#endif



