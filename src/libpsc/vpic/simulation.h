
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag
#include "util/rng/rng.h"

#include "field_array.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class Rng

#define IN_rng
#include "util/rng/rng_private.h"

struct Rng : rng {
  double drand();
  double drandn();
  double uniform(double lo, double hi);
  double normal(double mu, double sigma);
};

// ----------------------------------------------------------------------
// Rng implementation

inline double Rng::drand()
{
  return ::drand(this);
}

inline double Rng::drandn()
{
  return ::drandn(this);
}

inline double Rng::uniform(double lo, double hi)
{
  double dx = drand();
  return lo*(1.-dx) + hi*dx;
}

inline double Rng::normal(double mu, double sigma)
{
  return mu + sigma*drandn();
}

// ======================================================================
// class RngPool

struct RngPool {
  RngPool();
  
  void seed(int base, int which);
  Rng *operator[](int n);

  rng_pool *rng_pool_;
};

// ----------------------------------------------------------------------
// RngPool implementation

inline RngPool::RngPool()
{
  //rng_pool_ = simulation->entropy;
  int new_rng = 2;
  rng_pool_ = new_rng_pool(new_rng, 0, 0);
}

inline void RngPool::seed(int base, int which)
{
  seed_rng_pool(rng_pool_, base, which);
}

inline Rng* RngPool::operator[](int n)
{
  return reinterpret_cast<Rng *>(rng_pool_->rng[n]);
}

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
  
  RngPool rng_pool;

  //private:
  globals_diag *pDiag_;
  Grid grid_;
  MaterialList material_list_;
  field_array_t*& field_array_;
  interpolator_array_t*& interpolator_array_;
  accumulator_array_t*& accumulator_array_;
  hydro_array_t*& hydro_array_;
};

inline Simulation::Simulation()
  : grid_(simulation->grid),
    material_list_(simulation->material_list),
    field_array_(simulation->field_array),
    interpolator_array_(simulation->interpolator_array),
    accumulator_array_(simulation->accumulator_array),
    hydro_array_(simulation->hydro_array)
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


#endif



