
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag
#include "util/rng/rng.h"

// ======================================================================
// class Rng

#define IN_rng
#include "util/rng/rng_private.h"

struct Rng : rng {
  double uniform(double lo, double hi);
  double normal(double mu, double sigma);
};

// ----------------------------------------------------------------------
// Rng implementation

inline double Rng::uniform(double lo, double hi)
{
  return simulation->uniform(this, lo, hi);
}

inline double Rng::normal(double mu, double sigma)
{
  return simulation->normal(this, mu, sigma);
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

  RngPool rng_pool;
  //private:
  globals_diag *pDiag_;
};

// ----------------------------------------------------------------------
// Simulation implementation

inline Simulation::Simulation()
{
}

inline Simulation::~Simulation()
{
  delete pDiag_;
  // don't delete rng_pool_, because it's not ours
}

#endif



