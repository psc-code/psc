
#ifndef PSC_RNG_H
#define PSC_RNG_H

#include <random>
#include <cassert>
#include "mrc_common.h"
#include "psc_vpic_bits.h"

// ======================================================================
// PscRng

struct PscRng
{
  typedef std::mt19937 Urng;
  typedef std::uniform_real_distribution<double> Uniform;
  typedef std::normal_distribution<double> Normal;
  
  static PscRng* create()
  {
    return new PscRng;
  }

  void seed(unsigned int seed)        { urng_.seed(seed); }
  unsigned int operator()()           { return urng_(); }
  static constexpr unsigned int min() { return Urng::min(); }
  static constexpr unsigned int max() { return Urng::max(); }

  double uniform(double lo, double hi)
  {
    Uniform::param_type prm(lo, hi);
    return uniform_(urng_, prm);
  }
  
  double normal(double mu, double sigma)
  {
    Normal::param_type prm(mu, sigma);
    return normal_(urng_, prm);
  }
  
private:
  Urng urng_;
  Uniform uniform_;
  Normal normal_;
};

// ======================================================================
// PscRngPool

template<class R>
struct PscRngPool
{
  typedef R Rng;
  
  PscRngPool() :
    rng_(Rng::create()),
    n_rng_(2) // not really needed, but kept to keep seeds the same as vpic
  {
  }
  
  void seed(int base, int which)
  {
    assert(which == 0);
    int seed = psc_world_rank + (psc_world_size+1) * n_rng_ * base;
    rng_->seed(seed);
  }
  
  Rng *operator[](int n)
  {
    assert(n == 0);
    return rng_;
  }

  Rng *rng_;
  int n_rng_;
};

#endif
