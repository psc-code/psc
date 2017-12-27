
#ifndef VPIC_RNG_H
#define VPIC_RNG_H

#include <limits>
#include <cassert>

// ======================================================================
// Rng
//
// needs to support:
// static Rng* create()  (factory)
// void seed(unsigned int) 
// double uniform(double lo, double hi)
// double normal(double mu, double sigma)
//
// also supports UniformRandomBitGenerator concept, though that ends up
// only being tested but not used

// ======================================================================
// VpicRng

#define IN_rng
#include "util/rng/rng_private.h"

struct VpicRng : rng_t
{
  static VpicRng* create(int seed = 0)
  {
    return static_cast<VpicRng*>(::new_rng(seed));
  }
  
  void seed(unsigned int s)
  {
    ::seed_rng(this, s);
  }
  
  double uniform(double lo, double hi)
  {
    double dx = ::drand(this);
    return lo * (1.-dx) + hi * dx;
  }

  double normal(double mu, double sigma)
  {
    return mu + sigma * ::drandn(this);
  }

  unsigned int operator()()
  {
    return ::uirand(this);
  }

  static constexpr unsigned int min()
  {
    return std::numeric_limits<unsigned int>::min();
  }

  static constexpr unsigned int max()
  {
    return std::numeric_limits<unsigned int>::max();
  }
};

// ======================================================================
// RngPool
//
// needs to support:
// void seed(int base, int which)
// Rng* operator[](int n)

// ======================================================================
// VpicRngPool

template<class R>
struct VpicRngPool
{
  typedef R Rng;
  
  VpicRngPool()
  {
    int new_rng = 2;
    rng_pool_ = new_rng_pool(new_rng, 0, 0);
  }
  
  void seed(int base, int which)
  {
    ::seed_rng_pool(rng_pool_, base, which);
  }
  
  Rng* operator[](int n)
  {
    assert(n == 0);
    return reinterpret_cast<Rng *>(rng_pool_->rng[n]);
  }

  rng_pool *rng_pool_;
};


#endif


