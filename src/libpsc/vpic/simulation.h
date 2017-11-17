
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag
#include "util/rng/rng.h"

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
// Grid

struct Grid {
  Grid(grid *grid);
  
  void setup(double dx[3], double dt, double cvac, double eps0);
  void partition_periodic_box(double xl[3], double xh[3], int gdims[3], int np[3]);
  void set_fbc(int boundary, int fbc);
  void set_pbc(int boundary, int pbc);

private:
  grid *g_;
};

inline void Grid::setup(double dx[3], double dt, double cvac, double eps0)
{
  g_->dx = dx[0];
  g_->dy = dx[1];
  g_->dz = dx[2];
  g_->dt = dt;
  g_->cvac = cvac;
  g_->eps0 = eps0;
}

inline Grid::Grid(grid *g) :
  g_(g)
{
}

inline void Grid::partition_periodic_box(double xl[3], double xh[3],
					 int gdims[3], int np[3])
{
  ::partition_periodic_box(g_, xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
			   gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
}

inline void Grid::set_fbc(int boundary, int fbc)
{
  ::set_fbc(g_, boundary, fbc);
}

inline void Grid::set_pbc(int boundary, int pbc)
{
  ::set_pbc(g_, boundary, pbc);
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
  inline void define_field_array(struct field_array *fa, double damp);
  
  RngPool rng_pool;

  //private:
  globals_diag *pDiag_;
  Grid grid_;
};

// ----------------------------------------------------------------------
// Simulation implementation

inline Simulation::Simulation()
  : grid_(simulation->grid)
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
						    double eps, double mu,
						    double sigma, double zeta)
{
  return simulation->define_material(name, eps, mu, sigma, zeta);
}

inline void Simulation::define_field_array(struct field_array *fa, double damp)
{
  simulation->define_field_array(fa, damp);
}


#endif



