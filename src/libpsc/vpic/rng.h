
#ifndef RNG_H
#define RNG_H

// ======================================================================
// class Rng

#define IN_rng
#include "util/rng/rng_private.h"

struct VpicRng : rng {
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
};

// ======================================================================
// class RngPool

template<class Rng>
struct VpicRngPool
{
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
