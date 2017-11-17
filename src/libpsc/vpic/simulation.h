
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
// class Simulation

struct Simulation {
  Simulation();
  ~Simulation();

  void setup_grid(double dx[3], double dt, double cvac, double eps0);
  void define_periodic_grid(double xl[3], double xh[3], int gdims[3],
			    int np[3]);
  void set_domain_field_bc(struct Simulation *sim, int boundary, int bc);
  void set_domain_particle_bc(struct Simulation *sim, int boundary, int bc);
  
  RngPool rng_pool;
  //private:
  globals_diag *pDiag_;
  grid *grid_;
};

// ----------------------------------------------------------------------
// Simulation implementation

inline Simulation::Simulation()
{
  grid_ = simulation->grid;
}

inline Simulation::~Simulation()
{
  delete pDiag_;
}

inline void Simulation::setup_grid(double dx[3], double dt, double cvac, double eps0)
{
  grid_->dx = dx[0];
  grid_->dy = dx[1];
  grid_->dz = dx[2];
  grid_->dt = dt;
  grid_->cvac = cvac;
  grid_->eps0 = eps0;
}

inline void Simulation::define_periodic_grid(double xl[3], double xh[3], int gdims[3],
					     int np[3])
{
  simulation->px = size_t(np[0]);
  simulation->py = size_t(np[1]);
  simulation->pz = size_t(np[2]);
  partition_periodic_box(grid_, xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
			 gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
}

inline void Simulation::set_domain_field_bc(struct Simulation *sim, int boundary, int bc)
{
  int fbc;
  switch (bc) {
  case BND_FLD_CONDUCTING_WALL: fbc = pec_fields   ; break;
  case BND_FLD_ABSORBING:       fbc = absorb_fields; break;
  default: assert(0);
  }
  simulation->set_domain_field_bc(boundary, fbc);
}

inline void Simulation::set_domain_particle_bc(struct Simulation *sim, int boundary, int bc)
{
  int pbc;
  switch (bc) {
  case BND_PART_REFLECTING: pbc = reflect_particles; break;
  case BND_PART_ABSORBING:  pbc = absorb_particles ; break;
  default: assert(0);
  }
  simulation->set_domain_particle_bc(boundary, pbc);
}

#endif



