
#ifndef RNG_H
#define RNG_H

#include <limits>

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
  
  unsigned int uirand() { return ::uirand(this); }
  double drand() { return ::drand(this); }
  double drandn() { return ::drandn(this); }

  double uniform(double lo, double hi)
  {
    double dx = drand();
    return lo * (1.-dx) + hi * dx;
  }

  double normal(double mu, double sigma)
  {
    return mu + sigma * drandn();
  }

  unsigned int operator()()
  {
    return uirand();
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
  
  Rng *operator[](int n)
  {
    return reinterpret_cast<Rng *>(rng_pool_->rng[n]);
  }

  rng_pool *rng_pool_;
};

#endif
