
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag
#include "util/rng/rng.h"

// ----------------------------------------------------------------------
// class RngPool

struct RngPool {
  void seed(int base);
  Rng *operator[](int n);
};

// ----------------------------------------------------------------------
// class Simulation

struct Simulation {
  Simulation();
  ~Simulation();

  RngPool rng_pool;
  //private:
  globals_diag *pDiag_;
};

// ----------------------------------------------------------------------
// class Rng

#define IN_rng
#include "util/rng/rng_private.h"

struct Rng : rng {
  double uniform(double lo, double hi);
  double normal(double mu, double sigma);
};

// ======================================================================
// RngPool implementation

inline void RngPool::seed(int base)
{
  simulation->seed_entropy(base);
}

inline Rng* RngPool::operator[](int n)
{
  return static_cast<Rng *>(simulation->rng(n));
}

// ======================================================================
// Simulation implementation

inline Simulation::Simulation()
{
}

inline Simulation::~Simulation()
{
  delete pDiag_;
  // don't delete rng_pool_, because it's not ours
}

// ======================================================================
// Rng implementation

inline double Rng::uniform(double lo, double hi)
{
  return simulation->uniform(this, lo, hi);
}

inline double Rng::normal(double mu, double sigma)
{
  return simulation->normal(this, mu, sigma);
}

#endif



